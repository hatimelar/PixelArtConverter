importScripts(
  "webgl-utils.js",
  "bilateralShader.js",
  "pixelartShader.js",
  "intermediateFramebuffer.js"
);
onmessage = (event) => {
  console.log(event);
  if (event.data.job === "initialRender") {
    initialRender(event.data);
  } else if (event.data.job === "render") {
    createPixelart(event.data);
  } else if (event.data.job == "smoothAndRender") {
    smoothAndRenderPixelart(event.data);
  } else if (event.data.job === "startup") {
    startup();
    postMessage({ stage: "ready" });
  }
};

let offscreenCanvas;
let gl;
let image;
let smoothedImageData;
let bilateralShadingData;
let knollPixelartShadingData;
let orderedPixelartShadingData;
let intermediateBuffer;
let palettes;

function startup() {
  offscreenCanvas = new OffscreenCanvas(0, 0);
  gl = offscreenCanvas.getContext("webgl2");
  bilateralShadingData = setupBilateralFilterShadingProgram(gl);
  knollPixelartShadingData = setupPixelartShadingProgram(gl, true);
  orderedPixelartShadingData = setupPixelartShadingProgram(gl, false);
  intermediateBuffer = getIntermediateFramebuffer(gl);
}

function initialRender(msg) {
  postMessage({ stage: "Pre-processing image" });

  gl.canvas.width = msg.width;
  gl.canvas.height = msg.height;

  image = msg.image;

  palettes = generateColorPalettes(
    smoothImage(0, msg.image),
    msg.maxPaletteSize
  );

  smoothedImageData = smoothImage(msg.noiseReductionLevel, msg.image);

  createPixelart(msg);
}

function smoothAndRenderPixelart(msg) {
  postMessage({ stage: "Pre-processing image" });
  smoothedImageData = smoothImage(msg.noiseReductionLevel, image);
  createPixelart(msg);
}

function createPixelart(msg) {
  postMessage({ stage: "Rendering image" });
  const shadingData =
    msg.ditheringMethod === "Knoll"
      ? knollPixelartShadingData
      : orderedPixelartShadingData;
  const palette = msg.autoPalette
    ? palettes[msg.autoPaletteIndex]
    : msg.customPalette;

  renderPixelart(
    gl,
    shadingData,
    new PixelartConfig(
      msg.pixelSize,
      msg.ditheringMethod,
      msg.ditheringIntensity / 100,
      palette
    ),
    intermediateBuffer[1]
  );

  let imageData = getImageDataFromBuffer(null);

  postMessage({ stage: "Drawing edges" });
  /*if (msg.enableEdges) {
    let pixelMatrix = pixelMatrixFromImageData(smoothedImageData);
    pixelMatrix = subsamplePixelmatrix(pixelMatrix, msg.pixelSize);

    const edgeMatrix = cannyEdgeDetection(pixelMatrix, msg.edgeSensitivity);
    drawEdges(imageData, edgeMatrix, msg.pixelSize);
  }*/

  postMessage(imageData);
}

function smoothImage(iterations, image) {
  bilateralSmoothing(
    gl,
    bilateralShadingData,
    image,
    intermediateBuffer[0],
    intermediateBuffer[1],
    iterations
  );

  return getImageDataFromBuffer(intermediateBuffer[0]);
}

function getImageDataFromBuffer(fbo) {
  const pixels = new Uint8Array(
    offscreenCanvas.width * offscreenCanvas.height * 4
  );

  gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
  gl.readPixels(
    0,
    0,
    offscreenCanvas.width,
    offscreenCanvas.height,
    gl.RGBA,
    gl.UNSIGNED_BYTE,
    pixels
  );
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);

  return new ImageData(
    new Uint8ClampedArray(pixels),
    offscreenCanvas.width,
    offscreenCanvas.height
  );
}

const K_X = [
  [-1, 0, 1],
  [-2, 0, 2],
  [-1, 0, 1],
];

const K_Y = [
  [-1, -2, -1],
  [0, 0, 0],
  [1, 2, 1],
];

class Color {
  constructor(r, g, b, a) {
    this.r = r;
    this.g = g;
    this.b = b;
    this.a = a;
  }
}
const COLORS = {
  RED: "RED",
  GREEN: "GREEN",
  BLUE: "BLUE",
};

class Interval {
  constructor(low, high) {
    this.low = low;
    this.high = high;
  }
}

function clamp(number, low, high) {
  return Math.max(Math.min(number, high), low);
}

function generateColorPalettes(imageData, maxPaletteSize) {
  postMessage({ stage: "Scanning image" });
  let colors = getUniqueColors(imageData);
  postMessage({ stage: "Creating color palletes" });

  let palettes = medianCut(colors, maxPaletteSize);
  postMessage(palettes);
  return palettes;
}

function convertDecToTwoDigitHex(num) {
  let res = num.toString(16);
  return res.length == 1 ? "0" + res : res;
}

function convertColorArrToHexArr(colors) {
  const hexArr = [];
  for (const color of colors) {
    hexStr = `#${convertDecToTwoDigitHex(color.r)}${convertDecToTwoDigitHex(
      color.g
    )}${convertDecToTwoDigitHex(color.b)}`;
    hexArr.push(hexStr);
  }
  return hexArr;
}

function cannyEdgeDetection(pixelMatrix, edgeSensitivity) {
  pixelMatrix = getGrayscalePixelmatrix(pixelMatrix);

  let { gradient: g, theta: theta } = getGradientAndTheta(pixelMatrix);

  g = nonMaximalSuppression(g, theta);

  const maxRatio =
    0.7 * (1 - edgeSensitivity / 100) + (0.05 * edgeSensitivity) / 100;
  const minRatio =
    0.2 * (1 - edgeSensitivity / 100) + (0.05 * edgeSensitivity) / 100;

  doubleThreshold(g, minRatio, maxRatio);
  hysteresis(g);

  return g;
}

function getGrayscalePixelmatrix(pixelMatrix) {
  const grayscalePixelMatrix = [];
  for (let y = 0; y < pixelMatrix.length; y++) {
    grayscalePixelMatrix.push([]);
    for (let x = 0; x < pixelMatrix[0].length; x++) {
      grayscalePixelMatrix[y][x] =
        (299 * pixelMatrix[y][x].r +
          587 * pixelMatrix[y][x].g +
          114 * pixelMatrix[y][x].b) /
        1000.0;
    }
  }
  return grayscalePixelMatrix;
}

function pixelMatrixFromImageData(imageData) {
  const pixelMatrix = Array.from({ length: imageData.height }, () =>
    Array(imageData.width)
  );

  for (let y = 0; y < imageData.height; y++) {
    for (let x = 0; x < imageData.width; x++) {
      const index = x * 4 + y * 4 * imageData.width;
      pixelMatrix[y][x] = new Color(
        imageData.data[index],
        imageData.data[index + 1],
        imageData.data[index + 2],
        imageData.data[index + 3]
      );
    }
  }
  return pixelMatrix;
}
function subsamplePixelmatrix(pixelMatrix, pixelSize) {
  let outImageBlocksX = Math.floor(pixelMatrix[0].length / pixelSize);
  let outImageBlocksY = Math.floor(pixelMatrix.length / pixelSize);

  const resultMatrix = Array.from({ length: outImageBlocksY }, () =>
    Array(outImageBlocksX)
  );

  for (let y = 0; y < outImageBlocksY; y++) {
    for (let x = 0; x < outImageBlocksX; x++) {
      resultMatrix[y][x] = pixelMatrix[y * pixelSize][x * pixelSize];
    }
  }
  return resultMatrix;
}

// Compares two colors in the following descending order of channels: red, green, blue
// Returns true if a < b and false otherwise
function colorComparator(a, b) {
  if (a.r != b.r) return a.r - b.r;
  if (a.g != b.g) return a.g - b.g;
  return a.b - b.b;
}

// Returns an array (Color[]) containing unique colors extracted from the image
function getUniqueColors(imageData) {
  const colors = [];
  const colorSet = new Set();

  const imgDataArr = imageData.data;

  // Adding pixels to a temporary vector
  let color = [];

  for (let i = 0; i < imgDataArr.length; i += 4) {
    color = new Color(
      imgDataArr[i],
      imgDataArr[i + 1],
      imgDataArr[i + 2],
      imgDataArr[i + 3]
    );

    if (!colorSet.has(`${color.r}${color.g}${color.b}${color.a}`)) {
      colorSet.add(`${color.r}${color.g}${color.b}${color.a}`);
      colors.push(color);
    }
  }
  return colors;
}

// Takes an array of colors and a begin and end index
// Returns an enum representing the color/channel with the highest range
function colorHighestRange(colors, begin, end) {
  // Finding the color with the highest range

  let maxR = colors[begin].r;
  let minR = colors[begin].r;

  let maxG = colors[begin].g;
  let minG = colors[begin].g;

  let maxB = colors[begin].b;
  let minB = colors[begin].b;
  let color = {};
  for (let i = begin; i < end; i++) {
    color = colors[i];
    if (color.r > maxR) maxR = color.r;
    if (color.r < minR) minR = color.r;

    if (color.g > maxG) maxG = color.g;
    if (color.g < minG) minG = color.g;

    if (color.b > maxB) maxB = color.b;
    if (color.b < minB) minB = color.b;
  }
  let rRange = maxR - minR;
  let gRange = maxG - minG;
  let bRange = maxB - minB;

  // If the red channels has the largest range
  if (rRange >= gRange && rRange >= bRange) {
    return COLORS.RED;
  }
  // If the green channel has the largest range
  else if (gRange >= rRange && gRange >= bRange) {
    return COLORS.GREEN;
  }
  // If the blue channel has the largest range
  else {
    return COLORS.BLUE;
  }
}
// Color[] colors, int begin, int end, int depth, Color[][] partitionedColors
// Takes an array of colors a beginning and end index and a two dimensional array of colors
// that gets populated with color partitions. The functions returns nothing.
const sortSegment = (array, low, high, comparator) => {
  const segment = array.slice(low, high);
  segment.sort(comparator);
  for (let i = low; i < high; i++) array[i] = segment[i - low];
};

const getIntervalSum = (colors, low, high) => {
  let colorSum = new Color(0, 0, 0, 1);

  for (let i = low; i < high; i++) {
    colorSum.r += colors[i].r;
    colorSum.g += colors[i].g;
    colorSum.b += colors[i].b;
  }
  return colorSum;
};

function medianCut(colors, maxPaletteSize) {
  let low = 0;
  let high = colors.length;
  let partitions = [new Interval(low, high)];

  let palettes = [];
  for (let itr = 2; itr <= maxPaletteSize; itr *= 2) {
    let newPartitions = [];

    for (let i = 0; i < partitions.length; i++) {
      let interval = partitions[i];

      let highestRange = colorHighestRange(colors, interval.low, interval.high);
      switch (highestRange) {
        case COLORS.RED:
          sortSegment(colors, interval.low, interval.high, rColorComparator);
          break;
        case COLORS.GREEN:
          sortSegment(colors, interval.low, interval.high, gColorComparator);
          break;
        case COLORS.BLUE:
          sortSegment(colors, interval.low, interval.high, bColorComparator);
      }

      let medianIndex =
        Math.floor((interval.high - interval.low) / 2) + interval.low;

      if (medianIndex == interval.low) {
        newPartitions = partitions;
        break;
      }

      newPartitions.push(new Interval(interval.low, medianIndex));
      newPartitions.push(new Interval(medianIndex, interval.high));
    }
    partitions = newPartitions;

    const partitionSum = [];
    for (let i = 0; i < partitions.length; i++) {
      partitionSum.push(
        getIntervalSum(colors, partitions[i].low, partitions[i].high)
      );
    }
    // This is really, really bad and unreadable but it works
    for (let i = itr <= 16 ? partitions.length / 2 - 1 : 0; i >= 0; i--) {
      let palette = [];
      let flag = 0;
      for (; flag < i && itr <= 16 && i != 0; flag++) {
        const color = new Color(0, 0, 0, 255);

        let index =
          flag % 2 == 0
            ? flag
            : partitionSum.length - 2 - Math.floor(flag / 2) * 2;
        color.r = Math.floor(
          (partitionSum[index].r + partitionSum[index + 1].r) /
            (partitions[index + 1].high - partitions[index].low + 1)
        );
        color.g = Math.floor(
          (partitionSum[index].g + partitionSum[index + 1].g) /
            (partitions[index + 1].high - partitions[index].low + 1)
        );
        color.b = Math.floor(
          (partitionSum[index].b + partitionSum[index + 1].b) /
            (partitions[index + 1].high - partitions[index].low + 1)
        );
        palette.push(color);
      }

      for (
        j = flag % 2 == 0 ? flag : flag + 1;
        j < partitions.length - Math.floor(flag / 2) * 2;
        j++
      ) {
        const color = new Color(0, 0, 0, 255);
        color.r = Math.floor(
          partitionSum[j].r / (partitions[j].high - partitions[j].low + 1)
        );
        color.g = Math.floor(
          partitionSum[j].g / (partitions[j].high - partitions[j].low + 1)
        );
        color.b = Math.floor(
          partitionSum[j].b / (partitions[j].high - partitions[j].low + 1)
        );
        palette.push(color);
      }

      palettes.push(convertColorArrToHexArr(palette));
    }
  }

  return palettes;
}

function rColorComparator(a, b) {
  return a.r - b.r;
}
function gColorComparator(a, b) {
  return a.g - b.g;
}
function bColorComparator(a, b) {
  return a.b - b.b;
}

// returns the maximum distance between two consecutive numbers in
// an array
function maxDistanceInArray(channel) {
  channel.sort((a, b) => {
    return a - b;
  });
  let distance = 0;
  let max = 0;
  for (let i = 1; i < channel.length; i++) {
    distance = Math.abs(channel[i] - channel[i - 1]);
    if (distance > max) max = distance;
  }
  return max;
}

function luminanceComparator(a, b) {
  const aLuminance =
    (299.0 * a.r + 587.0 * a.g + 114.0 * a.b) / (255.0 * 1000.0);
  const bLuminance =
    (299.0 * b.r + 587.0 * b.g + 114.0 * b.b) / (255.0 * 1000.0);

  return aLuminance - bLuminance;
}

// Assumes a square odd dimensional kernel
function convolve(arr, kernel) {
  const kernelCenter = (kernel.length - 1) / 2;
  const res = Array.from({ length: arr.length }, () =>
    Array(arr[0].length).fill(0)
  );

  for (let y = kernelCenter; y < arr.length - kernelCenter; y++) {
    for (let x = kernelCenter; x < arr[0].length - kernelCenter; x++) {
      let val = 0;

      for (let y1 = 0; y1 < kernel.length; y1++) {
        for (let x1 = 0; x1 < kernel[0].length; x1++) {
          const yOffset = y1 - kernelCenter;
          const xOffset = x1 - kernelCenter;

          val += kernel[y1][x1] * arr[y + yOffset][x + xOffset];
        }
      }
      res[y][x] = val;
    }
  }
  return res;
}

function drawEdges(imageData, edgeMatrix, pixelSize) {
  const width = imageData.width;
  for (let y = 0; y < edgeMatrix.length; y++) {
    for (let x = 0; x < edgeMatrix[0].length; x++) {
      if (edgeMatrix[y][x] != 255) continue;

      console.log(pixelSize);
      const pixelIndex = (x + y * width) * pixelSize * 4;
      for (let j = 0; j < pixelSize; j++) {
        for (let i = 0; i < pixelSize; i++) {
          const offset = (j * width + i) * 4;

          imageData.data[pixelIndex + offset] = 0;
          imageData.data[pixelIndex + offset + 1] = 0;
          imageData.data[pixelIndex + offset + 2] = 0;
        }
      }
    }
  }
}

function getGradientAndTheta(pixelArr) {
  const xConv = convolve(pixelArr, K_X);
  const yConv = convolve(pixelArr, K_Y);

  const g = Array.from({ length: pixelArr.length }, () =>
    Array(pixelArr[0].length).fill(0)
  );
  const theta = Array.from({ length: pixelArr.length }, () =>
    Array(pixelArr[0].length).fill(0)
  );

  let max = 0;

  for (let y = 0; y < pixelArr.length; y++) {
    for (let x = 0; x < pixelArr[0].length; x++) {
      const hyp = Math.sqrt(xConv[y][x] ** 2 + yConv[y][x] ** 2);
      if (hyp > max) max = hyp;

      g[y][x] = hyp;
      theta[y][x] = Math.atan2(yConv[y][x], xConv[y][x]);
    }
  }

  for (let y = 0; y < g.length; y++) {
    for (let x = 0; x < g[0].length; x++) {
      g[y][x] = (g[y][x] / max) * 255;
    }
  }
  return { gradient: g, theta: theta };
}

function nonMaximalSuppression(gradArr, angArr) {
  const res = Array.from({ length: gradArr.length }, () =>
    Array(gradArr[0].length).fill(0)
  );

  for (let y = 1; y < gradArr.length - 1; y++) {
    for (let x = 1; x < gradArr[0].length - 1; x++) {
      let deg = (angArr[y][x] * 180) / Math.PI;
      if (deg < 0) deg += 180;

      let e1 = 0;
      let e2 = 0;

      if ((deg >= 0 && deg < 21.5) || (deg <= 180 && deg >= 157.5)) {
        e1 = gradArr[y][x - 1];
        e2 = gradArr[y][x + 1];
      } else if (deg >= 21.5 && deg < 67.5) {
        e1 = gradArr[y - 1][x - 1];
        e2 = gradArr[y + 1][x + 1];
      } else if (deg >= 67.5 && deg < 112.5) {
        e1 = gradArr[y - 1][x];
        e2 = gradArr[y + 1][x];
      } else {
        e1 = gradArr[y + 1][x - 1];
        e2 = gradArr[y - 1][x + 1];
      }

      if (gradArr[y][x] > e1 && gradArr[y][x] > e2) res[y][x] = gradArr[y][x];
      else res[y][x] = 0;
    }
  }
  return res;
}

function hysteresis(edgeMatrix, weak = 25, strong = 255) {
  for (let y = 1; y < edgeMatrix.length - 1; y++) {
    for (let x = 1; x < edgeMatrix[0].length - 1; x++) {
      if (edgeMatrix[y][x] == weak) {
        let success = false;
        for (let y1 = -1; y1 <= 1 && !success; y1++) {
          for (let x1 = -1; x1 <= 1; x1++) {
            if (edgeMatrix[y + y1][x + x1] == strong) {
              success = true;
              break;
            }
          }
        }

        if (success) {
          edgeMatrix[y][x] = strong;
        } else {
          edgeMatrix[y][x] = 0;
        }
      }
    }
  }
}

function doubleThreshold(gradArr, minRatio = 0.1, maxRatio = 0.2) {
  let max = 0;
  for (let y = 0; y < gradArr.length; y++) {
    for (let x = 0; x < gradArr[0].length; x++) {
      if (gradArr[y][x] > max) max = gradArr[y][x];
    }
  }

  const maxThreshold = max * maxRatio;
  const minThreshold = maxThreshold * minRatio;

  for (let y = 0; y < gradArr.length; y++) {
    for (let x = 0; x < gradArr[0].length; x++) {
      if (gradArr[y][x] >= maxThreshold) gradArr[y][x] = 255;
      else if (gradArr[y][x] >= minThreshold) gradArr[y][x] = 25;
      else gradArr[y][x] = 0;
    }
  }
}

onmessage = (event) => {
  postMessage({ stage: "HELLLOO???" });
  const bitMap = event.data.bitMap;
  if (event.data.job === "generateColorPalettes") {
    postMessage({ stage: "HELLLOO2222" });
    const palettes = generateColorPalettes(bitMap, event.data.maxPaletteSize);
    //postMessage(palettes);
  } else if (event.data.job === "generatePixelArt") {
    debounceGeneration(event, imageData);
  }
};

const debounceGeneration = (event, imageData) => {
  const res = generatePixelartCustomPalette(
    imageData,
    event.data.palette,
    event.data.ditheringIntensity,
    event.data.ditheringMethod,
    event.data.pixelSize,
    event.data.noiseReductionLevel,
    event.data.enableEdges,
    event.data.edgeSensitivity
  );

  self.postMessage(res);
};

const BAYER8 = [
  [0, 32, 8, 40, 2, 34, 10, 42],
  [48, 16, 56, 24, 50, 18, 58, 26],
  [12, 44, 4, 36, 14, 46, 6, 38],
  [60, 28, 52, 20, 62, 30, 54, 22],
  [3, 35, 11, 43, 1, 33, 9, 41],
  [51, 19, 59, 27, 49, 17, 57, 25],
  [15, 47, 7, 39, 13, 45, 5, 37],
  [63, 31, 55, 23, 61, 29, 53, 21],
];

const BAYER4 = [
  [0, 8, 2, 10],
  [12, 4, 14, 6],
  [3, 11, 1, 9],
  [15, 7, 13, 5],
];

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
function clamp(number, low, high) {
  return Math.max(Math.min(number, high), low);
}
function generateColorPalettes(imageBitmap, maxPaletteSize) {
  const offscreenCanvas = new OffscreenCanvas(
    imageBitmap.width,
    imageBitmap.height
  );
  const context = offscreenCanvas.getContext("2d");

  // Draw the ImageBitmap onto the OffscreenCanvas
  context.drawImage(imageBitmap, 0, 0);

  // Get ImageData from the OffscreenCanvas
  const imageData = context.getImageData(
    0,
    0,
    offscreenCanvas.width,
    offscreenCanvas.height
  );

  const palettes = [];
  postMessage({ stage: "Scanning image" });
  let colors = getUniqueColors(imageData);
  postMessage({ stage: "Creating color pallete" });
  postMessage({ stage: toString(colors.length) });
  for (let i = 2; i <= maxPaletteSize; i *= 2) {
    const partitionedColors = [];
    medianCut(colors, Math.log2(i), partitionedColors);
    const palette = getAvgColorFromSet(partitionedColors);
    palettes.push(convertColorArrToHexArr(palette));
  }
  postMessage({ stage: "Done" });
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

function generatePixelartAutoGenPalette(
  imageData,
  paletteSize,
  ditheringIntensity,
  pixelDimension
) {
  let partitionedColors = [];
  postMessage({ stage: "Scanning image" });
  let colors = getUniqueColors(imageData);
  postMessage({ stage: "Creating color pallete" });
  medianCut(colors, Math.log2(paletteSize), partitionedColors);
  let palette = getAvgColorFromSet(partitionedColors);
  postMessage({ stage: "Rendering image" });
  pixelizeOrdered(
    imageData,
    palette,
    8,
    pixelDimension,
    ditheringIntensity / 100
  );
}
function generatePixelartCustomPalette(
  imageData,
  hexArr,
  ditheringIntensity,
  ditheringMethod,
  pixelDimension,
  noiseReductionLevel,
  enableEdges,
  edgeSensitivity
) {
  postMessage({ stage: "Processing palette" });
  let palette = convertHexArrToPalette(hexArr);

  postMessage({ stage: "Rendering image" });
  const sig =
    (1 - noiseReductionLevel / 100) * 0 + (noiseReductionLevel / 100) * 4; // interpolate sigma between sig=1 and 4

  let pixelMatrix = pixelMatrixFromImageData(imageData);

  pixelMatrix = subsamplePixelmatrix(pixelMatrix, pixelDimension);

  if (ditheringMethod === "Knoll")
    pixelizeKnoll(pixelMatrix, palette, ditheringIntensity / 100);
  else pixelizeOrdered(pixelMatrix, palette, ditheringIntensity / 100, 8);

  if (enableEdges) {
    postMessage({ stage: "Detecting edges" });
    const edgeArr = getEdgeArr(imageData, pixelDimension, edgeSensitivity, sig);
    drawEdges(pixelMatrix, edgeArr);
  }
  return pixelArrToImageData(pixelMatrix, pixelDimension);
}
function getEdgeArr(imageData, pixelDimension, edgeSensitivity, sig = 2) {
  const pixelArr = getSmoothedGrayscale(imageData, pixelDimension, sig);

  let { gradient: g, theta: theta } = getGradientAndTheta(pixelArr);

  g = nonMaximalSuppression(g, theta);

  const maxRatio =
    0.7 * (1 - edgeSensitivity / 100) + (0.05 * edgeSensitivity) / 100;
  const minRatio =
    0.2 * (1 - edgeSensitivity / 100) + (0.05 * edgeSensitivity) / 100;
  doubleThreshold(g, minRatio, maxRatio);
  hysteresis(g);

  return g;
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

function bilateralSmoothing(
  pixelMatrix,
  diameter,
  center,
  g1Lookup,
  g2Lookup,
  vertical
) {
  const resultMatrix = Array.from({ length: pixelMatrix.length }, () =>
    Array(pixelMatrix[0].length)
  );

  for (let y = 0; y < pixelMatrix.length; y++) {
    for (let x = 0; x < pixelMatrix[0].length; x++) {
      let wR = 0;
      let wG = 0;
      let wB = 0;
      let filteredPixel = [0, 0, 0];

      if (vertical) {
        if (y < center || y > pixelMatrix.length - 1 - center) {
          resultMatrix[y][x] = pixelMatrix[y][x];
          continue;
        }
        for (let y1 = 0; y1 < diameter; y1++) {
          const j = y1 - center + y;
          const g1 = g1Lookup[Math.abs(y1 - center)];

          const cWR =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].r - pixelMatrix[j][x].r))
            ] * g1;
          const cWG =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].g - pixelMatrix[j][x].g))
            ] * g1;
          const cWB =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].b - pixelMatrix[j][x].b))
            ] * g1;

          wR += cWR;
          wG += cWG;
          wB += cWB;

          filteredPixel[0] += pixelMatrix[j][x].r * cWR;
          filteredPixel[1] += pixelMatrix[j][x].g * cWG;
          filteredPixel[2] += pixelMatrix[j][x].b * cWB;
        }
      } else {
        if (x < center || x > pixelMatrix[0].length - 1 - center) {
          resultMatrix[y][x] = pixelMatrix[y][x];
          continue;
        }
        for (let x1 = 0; x1 < diameter; x1++) {
          const i = x1 - center + x;
          const g1 = g1Lookup[Math.abs(x1 - center)];

          const cWR =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].r - pixelMatrix[y][i].r))
            ] * g1;
          const cWG =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].g - pixelMatrix[y][i].g))
            ] * g1;
          const cWB =
            g2Lookup[
              Math.round(Math.abs(pixelMatrix[y][x].b - pixelMatrix[y][i].b))
            ] * g1;

          wR += cWR;
          wG += cWG;
          wB += cWB;

          filteredPixel[0] += pixelMatrix[y][i].r * cWR;
          filteredPixel[1] += pixelMatrix[y][i].g * cWG;
          filteredPixel[2] += pixelMatrix[y][i].b * cWB;
        }
      }

      filteredPixel[0] /= wR;
      filteredPixel[1] /= wG;
      filteredPixel[2] /= wB;

      resultMatrix[y][x] = new Color(
        clamp(filteredPixel[0], 0, 255),
        clamp(filteredPixel[1], 0, 255),
        clamp(filteredPixel[2], 0, 255),
        255
      );
    }
  }
  return resultMatrix;
}

function createPaletteFromRGBA(rgbaArr) {
  let palette = [];
  for (let i = 0; i < rgbaArr.length; i += 4) {
    palette.push(
      new Color(rgbaArr[i], rgbaArr[i + 1], rgbaArr[i + 2], rgbaArr[i + 3])
    );
  }
  return palette;
}

function convertHexArrToPalette(hexArr) {
  const palette = [];

  for (let i = 0; i < hexArr.length; i++) {
    const hexStr = hexArr[i];

    const red = parseInt(hexStr[1], 16) * 16 + parseInt(hexStr[2], 16);
    const green = parseInt(hexStr[3], 16) * 16 + parseInt(hexStr[4], 16);
    const blue = parseInt(hexStr[5], 16) * 16 + parseInt(hexStr[6], 16);

    palette.push(new Color(red, green, blue, 255));
  }

  return palette;
}

// partitionedColors is of type Color[][]
// The function returns a pallete (Color[]) which corresponds
// to the union of the average color of each partition
function getAvgColorFromSet(partitionedColors) {
  // Linear color channel sums
  let palette = [];

  for (let partition of partitionedColors) {
    let rSum = 0;
    let gSum = 0;
    let bSum = 0;
    for (let color of partition) {
      rSum += color.r;
      gSum += color.g;
      bSum += color.b;
    }
    let partitionSize = partition.length;
    let r = Math.round(rSum / partitionSize);
    let g = Math.round(gSum / partitionSize);
    let b = Math.round(bSum / partitionSize);
    palette.push(new Color(r, g, b, 255));
  }
  return palette;
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
//Takes a color and an array of colors(pallete) and returns the color in the pallete that is the
// closest to it
function findNearestColor(mainColor, pallete) {
  let nearestColor;
  let minDistance = 999999;
  let distance;
  let rDiff;
  let gDiff;
  let bDiff;
  for (let color of pallete) {
    rDiff = color.r - mainColor.r;
    gDiff = color.g - mainColor.g;
    bDiff = color.b - mainColor.b;

    distance = 3 * (rDiff * rDiff) + 6 * (gDiff * gDiff) + bDiff * bDiff;
    if (distance < minDistance) {
      nearestColor = color;
      minDistance = distance;
    }
  }
  return nearestColor;
}

function findNearestColor2(mainColor, pallete) {
  let nearestColor;
  let minDistance = 999999;
  let distance;

  for (let color of pallete) {
    distance = Math.abs(luminanceComparator(mainColor, color));
    if (distance < minDistance) {
      nearestColor = color;
      minDistance = distance;
    }
  }
  return nearestColor;
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
  for (let i = begin; i <= end; i++) {
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
function medianCut(colors, depth, partitionedColors) {
  //Base case
  if (depth == 0 || colors.length == 1) {
    // Take the partition of pixels passed and add them to partitionedColors vector
    let partition = [];
    for (let i = 0; i < colors.length; i++) {
      partition.push(colors[i]);
    }
    partitionedColors.push(partition);
    return;
  }
  let highestRange = colorHighestRange(colors, 0, colors.length - 1);
  switch (highestRange) {
    case COLORS.RED:
      colors.sort(rColorComparator);
      break;
    case COLORS.GREEN:
      colors.sort(gColorComparator);
      break;
    case COLORS.BLUE:
      colors.sort(bColorComparator);
  }

  let medianIndex = Math.floor(colors.length / 2);

  medianCut(colors.slice(0, medianIndex), depth - 1, partitionedColors);
  medianCut(colors.slice(medianIndex), depth - 1, partitionedColors);
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
// Takes a pallete as input and returns an int[] representing the maximum
// difference between two consecutive values in each respective color channel
function maxDistanceChannel(pallete) {
  let result = [0, 0, 0];
  let rChannel = [];
  let gChannel = [];
  let bChannel = [];

  for (let i = 0; i < pallete.length; i++) {
    rChannel.push(pallete[i].r);
  }
  for (let i = 0; i < pallete.length; i++) {
    gChannel.push(pallete[i].g);
  }
  for (let i = 0; i < pallete.length; i++) {
    bChannel.push(pallete[i].b);
  }

  result[0] = maxDistanceInArray(rChannel);
  result[1] = maxDistanceInArray(gChannel);
  result[2] = maxDistanceInArray(bChannel);

  return result;
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

function pixelizeOrdered(pixelArr, pallete, ditherAmount, bayerSize) {
  let attempt;
  let outputColor;
  let r;
  let g;
  let b;

  const bayerScalarDenominator = bayerSize == 8 ? 64.0 : 16.0;
  let BAYER = [];
  if (bayerSize == 8) {
    BAYER = BAYER8;
    bayerSize = 8;
  } else {
    BAYER = BAYER4;
    bayerSize = 4;
  }

  const palleteOffset = maxDistanceChannel(pallete);

  for (let i = 0; i < pixelArr[0].length; i++) {
    for (let j = 0; j < pixelArr.length; j++) {
      r =
        pixelArr[j][i].r +
        palleteOffset[0] *
          ditherAmount *
          (BAYER[i % bayerSize][j % bayerSize] / bayerScalarDenominator - 0.5);
      g =
        pixelArr[j][i].g +
        palleteOffset[1] *
          ditherAmount *
          (BAYER[i % bayerSize][j % bayerSize] / bayerScalarDenominator - 0.5);
      b =
        pixelArr[j][i].b +
        palleteOffset[2] *
          ditherAmount *
          (BAYER[i % bayerSize][j % bayerSize] / bayerScalarDenominator - 0.5);

      attempt = new Color(
        clamp(Math.round(r), 0, 255),
        clamp(Math.round(g), 0, 255),
        clamp(Math.round(b), 0, 255),
        255
      );
      outputColor = findNearestColor(attempt, pallete);

      pixelArr[j][i] = { ...outputColor };
    }
  }
}
/*
function smoothImage(imageData, pixelSize, sig) {
  const mask = getGaussianMask(sig);
  const maskCenter = Math.floor((Math.ceil(sig * 3) * 2 - 1) / 2);

  const hPass = [];
  const res = [];

  let width = imageData.width;
  let outImageBlocksX = Math.floor(imageData.width / pixelSize);
  let outImageBlocksY = Math.floor(imageData.height / pixelSize);

  let imageDataArr = imageData.data;
  const pixelCenterOffsetX = Math.floor((pixelSize - 1) / 2) * 4;
  const pixelCenterOffsetY = pixelCenterOffsetX * width;

  for (let y = 0; y < outImageBlocksY; y++) {
    hPass.push([]);
    for (let x = 0; x < outImageBlocksX; x++) {
      inPosX = x * pixelSize;
      inPosY = y * pixelSize;

      // Pixel index is the center index of the superpixel
      let pixelIndex =
        (inPosX + inPosY * width) * 4 + pixelCenterOffsetX + pixelCenterOffsetY;

      const smoothedPixel = new Color(0, 0, 0, 255);

      let fail = false;
      for (let x1 = 0; x1 < mask.length; x1++) {
        const xOffset = (x1 - maskCenter) * 4;

        const tempPixelIndex = xOffset + pixelIndex;
        if (tempPixelIndex < 0 || tempPixelIndex >= imageDataArr.length) {
          fail = true;
          break;
        }

        smoothedPixel.r += imageDataArr[tempPixelIndex] * mask[x1];
        smoothedPixel.g += imageDataArr[tempPixelIndex + 1] * mask[x1];
        smoothedPixel.b += imageDataArr[tempPixelIndex + 2] * mask[x1];
      }

      if (!fail) {
        smoothedPixel.r = clamp(smoothedPixel.r, 0, 255);
        smoothedPixel.g = clamp(smoothedPixel.g, 0, 255);
        smoothedPixel.b = clamp(smoothedPixel.b, 0, 255);
      } else {
        smoothedPixel.r = imageDataArr[pixelIndex];
        smoothedPixel.g = imageDataArr[pixelIndex + 1];
        smoothedPixel.b = imageDataArr[pixelIndex + 2];
      }
      hPass[y].push(smoothedPixel);
    }
  }

  // Vertical pass
  for (let y = 0; y < hPass.length; y++) {
    res.push([]);
    for (let x = 0; x < hPass[0].length; x++) {
      let smoothedPixel = new Color(0, 0, 0, 255);

      for (let i = 0; i < mask.length; i++) {
        const yOffset = i - maskCenter;

        let y1 = yOffset + y;

        // Mirroring
        if (y1 < 0) {
          y1 = y1 * -1 - 1;
          const quo = Math.floor(y1 / hPass.length);

          y1 = y1 % hPass.length;

          if (quo % 2 != 0) {
            y1 = hPass.length - 1 - y1;
          }
        } else if (y1 >= hPass.length) {
          const quo = Math.floor(y1 / hPass.length);
          y1 = y1 % hPass.length;

          if (quo % 2 != 0) {
            y1 = hPass.length - 1 - y1;
          }
        }
        smoothedPixel.r += hPass[y1][x].r * mask[i];
        smoothedPixel.g += hPass[y1][x].g * mask[i];
        smoothedPixel.b += hPass[y1][x].b * mask[i];
      }

      res[y].push(smoothedPixel);
    }
  }

  return res;
}
*/

function luminanceComparator(a, b) {
  const aLuminance =
    (299.0 * a.r + 587.0 * a.g + 114.0 * a.b) / (255.0 * 1000.0);
  const bLuminance =
    (299.0 * b.r + 587.0 * b.g + 114.0 * b.b) / (255.0 * 1000.0);

  return aLuminance - bLuminance;
}

function pixelizeKnoll(pixelArr, palette, ditherAmount) {
  let potentialColors = [];
  let outputColor;

  //let channelMaxDist = maxDistanceChannel(palette);

  let index;

  for (let y = 0; y < pixelArr.length; y++) {
    for (let x = 0; x < pixelArr[0].length; x++) {
      const inputColor = pixelArr[y][x];
      const error = [0, 0, 0];
      potentialColors.length = 0;

      for (let i = 0; i < 16; i++) {
        let potentialColor = new Color(0, 0, 0, 255);

        potentialColor.r = clamp(
          inputColor.r + error[0] * ditherAmount,
          0,
          255
        );
        potentialColor.g = clamp(
          inputColor.g + error[1] * ditherAmount,
          0,
          255
        );
        potentialColor.b = clamp(
          inputColor.b + error[2] * ditherAmount,
          0,
          255
        );

        potentialColor = findNearestColor(potentialColor, palette);
        potentialColors.push(potentialColor);

        error[0] += inputColor.r - potentialColor.r;
        error[1] += inputColor.g - potentialColor.g;
        error[2] += inputColor.b - potentialColor.b;
        /*
        This fixes the underdithering when using custom palettes but can lead to more noise elsewhere
        if (
          error[0] * -1 >= channelMaxDist[0] * 4 &&
          (error[0] - inputColor.r + potentialColor.r) * -1 >=
            channelMaxDist[0] * 4
        )
          error[0] = channelMaxDist[0];
        else if (
          error[0] >= channelMaxDist[0] * 4 &&
          error[0] - inputColor.r + potentialColor.r >= channelMaxDist[0] * 4
        )
          error[0] = -channelMaxDist[0];

        if (
          error[1] * -1 >= channelMaxDist[1] * 4 &&
          (error[1] - inputColor.g + potentialColor.g) * -1 >=
            channelMaxDist[1] * 4
        )
          error[1] = channelMaxDist[1];
        else if (
          error[1] >= channelMaxDist[1] * 4 &&
          error[1] - inputColor.g + potentialColor.g >= channelMaxDist[1] * 4
        )
          error[1] = -channelMaxDist[1];

        if (
          error[2] * -1 >= channelMaxDist[2] * 4 &&
          (error[2] - inputColor.b + potentialColor.b) * -1 >=
            channelMaxDist[2] * 4
        )
          error[2] = channelMaxDist[2];
        else if (
          error[2] >= channelMaxDist[2] * 4 &&
          error[2] - inputColor.b + potentialColor.b >= channelMaxDist[2] * 4
        )
          error[2] = -channelMaxDist[2];
      */
      }

      potentialColors.sort(luminanceComparator);

      index = BAYER4[y % 4][x % 4];

      outputColor = potentialColors[index];

      pixelArr[y][x] = { ...outputColor };
    }
  }
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

function drawEdgesOnly(g, pixelSize) {
  const width = g[0].length * pixelSize;
  const height = g.length * pixelSize;
  const imageDataArr = new Uint8ClampedArray(width * height * 4);
  for (let i = 3; i < imageDataArr.length; i += 4) {
    imageDataArr[i] = 255;
  }

  for (let y = 0; y < g.length; y++) {
    for (let x = 0; x < g[0].length; x++) {
      const pixelOffset = x * pixelSize * 4 + y * pixelSize * 4 * width;

      for (let y1 = 0; y1 < pixelSize; y1++) {
        for (let x1 = 0; x1 < pixelSize; x1++) {
          const pixelIndex = pixelOffset + x1 * 4 + y1 * 4 * width;
          imageDataArr[pixelIndex] = g[y][x];
          imageDataArr[pixelIndex + 1] = g[y][x];
          imageDataArr[pixelIndex + 2] = g[y][x];
        }
      }
    }
  }
  return new ImageData(imageDataArr, width, height);
}

function drawEdges(pixelArr, edgeArr) {
  let c = 0;
  for (let y = 0; y < pixelArr.length; y++) {
    for (let x = 0; x < pixelArr[0].length; x++) {
      if (edgeArr[y][x] == 255) {
        pixelArr[y][x].r = 0;
        pixelArr[y][x].g = 0;
        pixelArr[y][x].b = 0;
        c++;
      }
    }
  }
}

function pixelArrToImageData(pixelArr, pixelSize) {
  const width = pixelArr[0].length * pixelSize;
  const height = pixelArr.length * pixelSize;
  const imageDataArr = new Uint8ClampedArray(width * height * 4);
  for (let i = 3; i < imageDataArr.length; i += 4) {
    imageDataArr[i] = 255;
  }

  for (let y = 0; y < pixelArr.length; y++) {
    for (let x = 0; x < pixelArr[0].length; x++) {
      const pixelOffset = x * pixelSize * 4 + y * pixelSize * 4 * width;

      for (let y1 = 0; y1 < pixelSize; y1++) {
        for (let x1 = 0; x1 < pixelSize; x1++) {
          const pixelIndex = pixelOffset + x1 * 4 + y1 * 4 * width;
          imageDataArr[pixelIndex] = pixelArr[y][x].r;
          imageDataArr[pixelIndex + 1] = pixelArr[y][x].g;
          imageDataArr[pixelIndex + 2] = pixelArr[y][x].b;
        }
      }
    }
  }
  return new ImageData(imageDataArr, width, height);
}

function getSmoothedGrayscale(imageData, pixelSize, sig) {
  const mask = sig != 0 ? getGaussianMask(sig) : [[1]];
  const maskCenter = (Math.ceil(sig * 2) * 2 + 1 - 1) / 2;

  const res = [];
  let width = imageData.width;
  let outImageBlocksX = Math.floor(imageData.width / pixelSize);
  let outImageBlocksY = Math.floor(imageData.height / pixelSize);

  let imageDataArr = imageData.data;
  const pixelCenterOffsetX = Math.floor((pixelSize - 1) / 2) * 4;
  const pixelCenterOffsetY = pixelCenterOffsetX * width;

  for (let y = 0; y < outImageBlocksY; y++) {
    res.push([]);
    for (let x = 0; x < outImageBlocksX; x++) {
      inPosX = x * pixelSize;
      inPosY = y * pixelSize;

      // Pixel index is the center index of the superpixel
      let pixelIndex =
        (inPosX + inPosY * width) * 4 + pixelCenterOffsetX + pixelCenterOffsetY;

      let smoothedPixel = 0;

      let fail = false;
      for (let y1 = 0; y1 < mask.length; y1++) {
        for (let x1 = 0; x1 < mask[0].length; x1++) {
          const xOffset = (x1 - maskCenter) * 4;
          const yOffset = (y1 - maskCenter) * width * 4;

          const tempPixelIndex = xOffset + yOffset + pixelIndex;
          if (tempPixelIndex < 0 || tempPixelIndex >= imageDataArr.length) {
            fail = true;
            continue;
          }
          smoothedPixel +=
            (imageDataArr[tempPixelIndex] * 0.299 +
              imageDataArr[tempPixelIndex + 1] * 0.587 +
              imageDataArr[tempPixelIndex + 2] * 0.114) *
            mask[y1][x1];
        }
      }

      if (!fail) {
        smoothedPixel = clamp(smoothedPixel, 0, 255);
      } else {
        smoothedPixel =
          imageDataArr[pixelIndex] * 0.299 +
          imageDataArr[pixelIndex + 1] * 0.587 +
          imageDataArr[pixelIndex + 2] * 0.114;
      }
      res[y].push(smoothedPixel);
    }
  }
  return res;
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

function smoothImage(imageData, pixelSize, sig) {
  const mask = sig != 0 ? getGaussianMask(sig) : [[1]];
  const maskCenter = (Math.ceil(sig * 2) * 2 + 1 - 1) / 2;

  const res = [];
  let width = imageData.width;
  let outImageBlocksX = Math.floor(imageData.width / pixelSize);
  let outImageBlocksY = Math.floor(imageData.height / pixelSize);

  let imageDataArr = imageData.data;
  const pixelCenterOffsetX = Math.floor((pixelSize - 1) / 2) * 4;
  const pixelCenterOffsetY = pixelCenterOffsetX * width;

  for (let y = 0; y < outImageBlocksY; y++) {
    res.push([]);
    for (let x = 0; x < outImageBlocksX; x++) {
      inPosX = x * pixelSize;
      inPosY = y * pixelSize;

      // Pixel index is the center index of the superpixel
      let pixelIndex =
        (inPosX + inPosY * width) * 4 + pixelCenterOffsetX + pixelCenterOffsetY;

      const smoothedPixel = new Color(0, 0, 0, 255);
      let fail = false;
      for (let y1 = 0; y1 < mask.length; y1++) {
        for (let x1 = 0; x1 < mask[0].length; x1++) {
          const xOffset = (x1 - maskCenter) * 4;
          const yOffset = (y1 - maskCenter) * width * 4;

          const tempPixelIndex = xOffset + yOffset + pixelIndex;
          if (tempPixelIndex < 0 || tempPixelIndex >= imageDataArr.length) {
            fail = true;
            continue;
          }

          smoothedPixel.r += imageDataArr[tempPixelIndex] * mask[y1][x1];
          smoothedPixel.g += imageDataArr[tempPixelIndex + 1] * mask[y1][x1];
          smoothedPixel.b += imageDataArr[tempPixelIndex + 2] * mask[y1][x1];
        }
      }

      if (!fail) {
        smoothedPixel.r = clamp(smoothedPixel.r, 0, 255);
        smoothedPixel.g = clamp(smoothedPixel.g, 0, 255);
        smoothedPixel.b = clamp(smoothedPixel.b, 0, 255);
      } else {
        smoothedPixel.r = imageDataArr[pixelIndex];
        smoothedPixel.g = imageDataArr[pixelIndex + 1];
        smoothedPixel.b = imageDataArr[pixelIndex + 2];
      }
      res[y].push(smoothedPixel);
    }
  }
  return res;
}

function boxBlur(imageData, pixelSize, dim) {
  const res = [];
  let width = imageData.width;
  let outImageBlocksX = Math.floor(imageData.width / pixelSize);
  let outImageBlocksY = Math.floor(imageData.height / pixelSize);

  let imageDataArr = imageData.data;

  for (let y = 0; y < outImageBlocksY; y++) {
    res.push([]);
    for (let x = 0; x < outImageBlocksX; x++) {
      inPosX = x * pixelSize;
      inPosY = y * pixelSize;

      // Pixel index is the center index of the superpixel
      let pixelIndex = (inPosX + inPosY * width) * 4;
      const smoothedPixel = new Color(0, 0, 0, 255);

      let fail = false;
      for (let y1 = 0; y1 < dim; y1++) {
        for (let x1 = 0; x1 < dim; x1++) {
          const xOffset = x1 * 4;
          const yOffset = y1 * width * 4;

          const tempPixelIndex = xOffset + yOffset + pixelIndex;
          if (tempPixelIndex < 0 || tempPixelIndex >= imageDataArr.length) {
            fail = true;
            continue;
          }

          smoothedPixel.r += imageDataArr[tempPixelIndex];
          smoothedPixel.g += imageDataArr[tempPixelIndex + 1];
          smoothedPixel.b += imageDataArr[tempPixelIndex + 2];
        }
      }

      if (!fail) {
        smoothedPixel.r = clamp(smoothedPixel.r / dim ** 2, 0, 255);
        smoothedPixel.g = clamp(smoothedPixel.g / dim ** 2, 0, 255);
        smoothedPixel.b = clamp(smoothedPixel.b / dim ** 2, 0, 255);
      } else {
        smoothedPixel.r = imageDataArr[pixelIndex];
        smoothedPixel.g = imageDataArr[pixelIndex + 1];
        smoothedPixel.b = imageDataArr[pixelIndex + 2];
      }
      res[y].push(smoothedPixel);
    }
  }
  return res;
}

// size of kernel/mask is Math.ceil(sig)*4-1 x Math.ceil(sig)*4-1
function getGaussianMask(sig) {
  const numPoints = Math.ceil(sig * 2) * 2 + 1;
  const smoothingMask = [];
  const center = (numPoints - 1) / 2;

  let sum = 0;
  for (let y = 0; y < numPoints; y++) {
    smoothingMask.push([]);
    for (let x = 0; x < numPoints; x++) {
      smoothingMask[y].push(gaussValAt(x, y, center, center, sig));
      sum += smoothingMask[y][x];
    }
  }

  // Normalize mask so it sums to 1
  for (let y = 0; y < numPoints; y++) {
    for (let x = 0; x < numPoints; x++) {
      smoothingMask[y][x] = smoothingMask[y][x] / sum;
    }
  }
  return smoothingMask;
}

// Evaluates a gaussian function specified by the given parameters
function gaussValAt(x, y, ux, uy, sig) {
  const exponentX = -((x - ux) ** 2) / (2 * sig ** 2);
  const exponentY = -((y - uy) ** 2) / (2 * sig ** 2);

  return (1 / (sig * Math.sqrt(2 * Math.PI))) * Math.exp(exponentX + exponentY);
}
function gaussValAt2(x, sig) {
  const exponentX = -(x ** 2) / (2 * sig ** 2);

  return (1 / (sig ** 2 * 2 * Math.PI)) * Math.exp(exponentX);
}

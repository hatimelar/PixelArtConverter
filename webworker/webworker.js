importScripts(
  "Color.js",
  "colorQuantizer.js",
  "edgeDetector.js",
  "webgl-utils.js",
  "bilateralShader.js",
  "pixelartShader.js",
  "intermediateFramebuffer.js"
);
onmessage = (event) => {
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
  gl.finish();
  let imageData = getImageDataFromBuffer(null);

  if (msg.enableEdges) {
    postMessage({ stage: "Drawing edges" });
    let pixelMatrix = pixelMatrixFromImageData(smoothedImageData);
    pixelMatrix = subsamplePixelmatrix(pixelMatrix, msg.pixelSize);

    const edgeMatrix = cannyEdgeDetection(pixelMatrix, msg.edgeSensitivity);
    drawEdges(imageData, edgeMatrix, msg.pixelSize);
  }

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
  gl.finish();
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

function generateColorPalettes(imageData, maxPaletteSize) {
  postMessage({ stage: "Scanning image" });
  let colors = getUniqueColors(imageData);
  postMessage({ stage: "Creating color palletes" });

  let palettes = medianCut(colors, maxPaletteSize);
  postMessage(palettes);
  return palettes;
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

function drawEdges(imageData, edgeMatrix, pixelSize) {
  const width = imageData.width;
  for (let y = 0; y < edgeMatrix.length; y++) {
    for (let x = 0; x < edgeMatrix[0].length; x++) {
      if (edgeMatrix[y][x] != 255) continue;

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

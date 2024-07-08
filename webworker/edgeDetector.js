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

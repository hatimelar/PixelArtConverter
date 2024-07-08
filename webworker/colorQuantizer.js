class Interval {
  constructor(low, high) {
    this.low = low;
    this.high = high;
  }
}

// Compares two colors in the following descending order of channels: red, green, blue
// Returns true if a < b and false otherwise
function colorComparator(a, b) {
  if (a.r != b.r) return a.r - b.r;
  if (a.g != b.g) return a.g - b.g;
  return a.b - b.b;
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

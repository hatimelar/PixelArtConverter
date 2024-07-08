"use strict";
var vertexShaderSource = `#version 300 es

// an attribute is an input (in) to a vertex shader.
// It will receive data from a buffer
in vec2 a_position;
in vec2 a_texCoord;

// Used to pass in the resolution of the canvas
uniform vec2 u_resolution;


// Used to pass the texture coordinates to the fragment shader
out vec2 v_texCoord;

// all shaders have a main function
void main() {

  // convert the position from pixels to 0.0 to 1.0
  vec2 zeroToOne = a_position / u_resolution;

  // convert from 0->1 to 0->2
  vec2 zeroToTwo = zeroToOne * 2.0;

  // convert from 0->2 to -1->+1 (clipspace)
  vec2 clipSpace = zeroToTwo - 1.0;

  gl_Position = vec4(clipSpace, 0, 1);

  // pass the texCoord to the fragment shader
  // The GPU will interpolate this value between points.
  v_texCoord = a_texCoord;
}
`;

var knollFragmentShaderSource = `#version 300 es
precision highp float;

uniform sampler2D u_image;
uniform float u_pixelSize;
uniform vec4 u_palette[256];
uniform int u_paletteSize;
uniform float u_ditheringIntensity;
uniform vec2 u_resolution;

uniform int u_BAYER4[16];

in vec2 v_texCoord;
out vec4 outColor;


float luminosityComparator(vec4 a, vec4 b){
  float aLuminance = 0.299 * a.r + 0.587 * a.g + 0.114 * a.b;
  float bLuminance = 0.299 * b.r + 0.587 * b.g + 0.114 * b.b;

  return aLuminance - bLuminance;
}

vec4 findNearestColor(vec4 color){
  int index =0;
  float minDistance = 1.0e6;

  float rDiff=0.;
  float gDiff=0.;
  float bDiff=0.;
  float diff =0.;


  for (int i =0;i<u_paletteSize;i++){
    rDiff = color.r - u_palette[i].r;
    gDiff = color.g - u_palette[i].g;
    bDiff = color.b - u_palette[i].b;

    diff = 3. * rDiff*rDiff + 6. * gDiff*gDiff + bDiff*bDiff;
    
    if (diff<minDistance){
      minDistance = diff;
      index = i;
    }    
  }
  return u_palette[index];
}

void main() {
  vec2 pixelCoord = v_texCoord* u_resolution;
  pixelCoord = u_pixelSize*floor(pixelCoord/u_pixelSize)/u_resolution;
  
  
  vec3 error = vec3(0.,0.,0.);
  vec4 potentialColors[16];
  vec4 initialColor = texture(u_image, pixelCoord);

  
  for (int i=0; i< 16; i++){
    vec4 potentialColor = vec4(0.,0.,0.,1.);
    potentialColor.r = clamp(initialColor.r +error[0]*u_ditheringIntensity, 0., 1.);
    potentialColor.g = clamp(initialColor.g +error[1]*u_ditheringIntensity, 0., 1.);
    potentialColor.b = clamp(initialColor.b +error[2]*u_ditheringIntensity, 0., 1.);

    potentialColor = findNearestColor(potentialColor);
    potentialColors[i] = potentialColor;

    error[0] += initialColor.r - potentialColor.r;
    error[1] += initialColor.g - potentialColor.g;
    error[2] += initialColor.b - potentialColor.b;
  }

  for (int i=0;i<16;i++){
    vec4 color = potentialColors[i];

    int j;
    for (j=i-1; j>=0 && luminosityComparator(potentialColors[j],color)>0.0; j--){
      potentialColors[j+1] = potentialColors[j];
    }
    potentialColors[j+1] = color;
  }
  ivec2 superPixelCoord = ivec2(v_texCoord * u_resolution/u_pixelSize);
  outColor = potentialColors[u_BAYER4[(superPixelCoord.x%4) +((superPixelCoord.y%4)*4)]];
}`;

var orderedFragmentShaderSource = `#version 300 es
precision highp float;

uniform sampler2D u_image;
uniform float u_pixelSize;
uniform vec4 u_palette[256];
uniform int u_paletteSize;
uniform float u_ditheringIntensity;
uniform vec2 u_resolution;

uniform int u_BAYER4[16];
uniform vec3 u_channelOffset;

in vec2 v_texCoord;
out vec4 outColor;

vec4 findNearestColor(vec4 color){
  int index =0;
  float minDistance = 1.0e6;

  float rDiff=0.;
  float gDiff=0.;
  float bDiff=0.;
  float diff =0.;


  for (int i =0;i<u_paletteSize;i++){
    rDiff = color.r - u_palette[i].r;
    gDiff = color.g - u_palette[i].g;
    bDiff = color.b - u_palette[i].b;

    diff = 3. * rDiff*rDiff + 6. * gDiff*gDiff + bDiff*bDiff;
    
    if (diff<minDistance){
      minDistance = diff;
      index = i;
    }    
  }
  return u_palette[index];
}

void main() {
  vec2 pixelCoord = v_texCoord* u_resolution;
  pixelCoord = u_pixelSize*floor(pixelCoord/u_pixelSize)/u_resolution;

  ivec2 superPixelCoord = ivec2(v_texCoord * u_resolution/u_pixelSize);

  vec4 initialColor = texture(u_image, pixelCoord);
  
  float bayerCoeffecient  = (float(u_BAYER4[superPixelCoord.x % 4 + ((superPixelCoord.y%4)*4)]) / 16.0) - 0.5;

  float r =
    initialColor.r +
    u_channelOffset.r *
    u_ditheringIntensity *
    bayerCoeffecient;

  float g =
    initialColor.g +
    u_channelOffset.g *
    u_ditheringIntensity *
    bayerCoeffecient;

  float b =
    initialColor.b +
    u_channelOffset.b *
    u_ditheringIntensity *
    bayerCoeffecient;
  
  outColor = findNearestColor(vec4(clamp(r,0.,1.),clamp(g,0.,1.),clamp(b,0.,1.), 1.));
          
}`;

class PixelartConfig {
  constructor(pixelSize, ditheringMethod, ditheringIntensity, palette) {
    this.pixelSize = pixelSize;
    this.ditheringMethod = ditheringMethod;
    this.ditheringIntensity = ditheringIntensity;

    this.palette = convertHexArrToRGBAArr(palette);
  }
}

class PixelartShadingData {
  constructor(
    program,
    positionAttributeLocation,
    texCoordAttributeLocation,
    resolutionLocation,
    imageLocation,
    colorPaletteLocation,
    paletteSizeLocation,
    pixelSizeLocation,
    ditheringIntensityLocation,
    bayerArrayLocation,
    vao,
    positionBuffer,
    originalImageTexture
  ) {
    this.program = program;
    this.positionAttributeLocation = positionAttributeLocation;
    this.texCoordAttributeLocation = texCoordAttributeLocation;
    this.resolutionLocation = resolutionLocation;
    this.imageLocation = imageLocation;
    this.vao = vao;
    this.positionBuffer = positionBuffer;
    this.originalImageTexture = originalImageTexture;
    this.colorPaletteLocation = colorPaletteLocation;
    this.paletteSizeLocation = paletteSizeLocation;
    this.pixelSizeLocation = pixelSizeLocation;
    this.ditheringIntensityLocation = ditheringIntensityLocation;
    this.bayerArrayLocation = bayerArrayLocation;
  }
}

function setupPixelartShadingProgram(gl, useKnollDithering) {
  // setup GLSL program
  var program = createProgramFromSources(gl, [
    vertexShaderSource,
    useKnollDithering ? knollFragmentShaderSource : orderedFragmentShaderSource,
  ]);

  // look up where the vertex data needs to go.
  var positionAttributeLocation = gl.getAttribLocation(program, "a_position");
  var texCoordAttributeLocation = gl.getAttribLocation(program, "a_texCoord");

  // lookup uniforms
  var resolutionLocation = gl.getUniformLocation(program, "u_resolution");
  var imageLocation = gl.getUniformLocation(program, "u_image");
  var colorPaletteLocation = gl.getUniformLocation(program, "u_palette");
  var paletteSizeLocation = gl.getUniformLocation(program, "u_paletteSize");
  var pixelSizeLocation = gl.getUniformLocation(program, "u_pixelSize");
  var ditheringIntensityLocation = gl.getUniformLocation(
    program,
    "u_ditheringIntensity"
  );
  var bayerArrayLocation = gl.getUniformLocation(program, "u_BAYER4");

  // Create a vertex array object (attribute state)
  var vao = gl.createVertexArray();

  // and make it the one we're currently working with
  gl.bindVertexArray(vao);

  // Create a buffer and put a single pixel space rectangle in
  // it (2 triangles)
  var positionBuffer = gl.createBuffer();

  // Turn on the attribute
  gl.enableVertexAttribArray(positionAttributeLocation);

  // Bind it to ARRAY_BUFFER (think of it as ARRAY_BUFFER = positionBuffer)
  gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);

  // Tell the attribute how to get data out of positionBuffer (ARRAY_BUFFER)
  var size = 2; // 2 components per iteration
  var type = gl.FLOAT; // the data is 32bit floats
  var normalize = false; // don't normalize the data
  var stride = 0; // 0 = move forward size * sizeof(type) each iteration to get the next position
  var offset = 0; // start at the beginning of the buffer
  gl.vertexAttribPointer(
    positionAttributeLocation,
    size,
    type,
    normalize,
    stride,
    offset
  );

  // provide texture coordinates for the rectangle.
  var texCoordBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
  gl.bufferData(
    gl.ARRAY_BUFFER,
    new Float32Array([
      0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0,
    ]),
    gl.STATIC_DRAW
  );

  // Turn on the attribute
  gl.enableVertexAttribArray(texCoordAttributeLocation);

  // Tell the attribute how to get data out of texCoordBuffer (ARRAY_BUFFER)
  var size = 2; // 2 components per iteration
  var type = gl.FLOAT; // the data is 32bit floats
  var normalize = false; // don't normalize the data
  var stride = 0; // 0 = move forward size * sizeof(type) each iteration to get the next position
  var offset = 0; // start at the beginning of the buffer
  gl.vertexAttribPointer(
    texCoordAttributeLocation,
    size,
    type,
    normalize,
    stride,
    offset
  );

  function createAndSetupTexture(gl) {
    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);

    // Set up texture so we can render any size image and so we are
    // working with pixels.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    return texture;
  }

  // Create a texture and put the image in it.
  var originalImageTexture = createAndSetupTexture(gl);

  return new PixelartShadingData(
    program,
    positionAttributeLocation,
    texCoordAttributeLocation,
    resolutionLocation,
    imageLocation,
    colorPaletteLocation,
    paletteSizeLocation,
    pixelSizeLocation,
    ditheringIntensityLocation,
    bayerArrayLocation,
    vao,
    positionBuffer,
    originalImageTexture
  );
}

function renderPixelart(gl, pixelartShadingData, pixelartConfig, image) {
  // Upload the image into the texture.
  gl.bindTexture(gl.TEXTURE_2D, pixelartShadingData.originalImageTexture);
  // Upload the image into the texture.
  var mipLevel = 0; // the largest mip
  var internalFormat = gl.RGBA; // format we want in the texture
  var srcFormat = gl.RGBA; // format of data we are supplying
  var srcType = gl.UNSIGNED_BYTE; // type of data we are supplying
  gl.texImage2D(
    gl.TEXTURE_2D,
    mipLevel,
    internalFormat,
    srcFormat,
    srcType,
    image
  );
  // Bind the position buffer so gl.bufferData that will be called
  // in setRectangle puts data in the position buffer
  gl.bindBuffer(gl.ARRAY_BUFFER, pixelartShadingData.positionBuffer);

  // Set a rectangle the same size as the image.
  setRectangle(gl, 0, 0, gl.canvas.width, gl.canvas.height);

  // Tell WebGL how to convert from clip space to pixels
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

  // Clear the canvas
  gl.clearColor(0, 0, 0, 0);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  // Tell it to use our program (pair of shaders)
  gl.useProgram(pixelartShadingData.program);

  // Bind the attribute/buffer set we want.
  gl.bindVertexArray(pixelartShadingData.vao);

  // start with the original image on unit 0
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, pixelartShadingData.originalImageTexture);

  // Tell the shader to get the texture from texture unit 0
  gl.uniform1i(pixelartShadingData.imageLocation, 0);

  // Set other uniforms
  gl.uniform4fv(
    pixelartShadingData.colorPaletteLocation,
    pixelartConfig.palette
  );
  gl.uniform1i(
    pixelartShadingData.paletteSizeLocation,
    pixelartConfig.palette.length / 4
  );
  gl.uniform1f(
    pixelartShadingData.pixelSizeLocation,
    parseFloat(pixelartConfig.pixelSize)
  );
  gl.uniform1f(
    pixelartShadingData.ditheringIntensityLocation,
    pixelartConfig.ditheringIntensity
  );

  // Flatened 4x4 Bayer matrix
  gl.uniform1iv(
    pixelartShadingData.bayerArrayLocation,
    new Int32Array([0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5])
  );

  if (pixelartConfig.ditheringMethod === "Ordered") {
    const channelOffsetLocation = gl.getUniformLocation(
      pixelartShadingData.program,
      "u_channelOffset"
    );

    gl.uniform3fv(
      channelOffsetLocation,
      maxDistanceChannel(pixelartConfig.palette)
    );
  }

  setFramebuffer(null, gl.canvas.width, gl.canvas.height);

  // Clear the canvas
  gl.clearColor(0, 0, 0, 0);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  draw();
  gl.bindTexture(gl.TEXTURE_2D, null);

  function draw() {
    // Draw the rectangle.
    var primitiveType = gl.TRIANGLES;
    var offset = 0;
    var count = 6;
    gl.drawArrays(primitiveType, offset, count);
  }

  function setRectangle(gl, x, y, width, height) {
    var x1 = x;
    var x2 = x + width;
    var y1 = y;
    var y2 = y + height;
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([x1, y1, x2, y1, x1, y2, x1, y2, x2, y1, x2, y2]),
      gl.STATIC_DRAW
    );
  }

  function setFramebuffer(fbo, width, height) {
    // make this the framebuffer we are rendering to.
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);

    // Tell the shader the resolution of the framebuffer.
    gl.uniform2f(pixelartShadingData.resolutionLocation, width, height);

    // Tell WebGL how to convert from clip space to pixels
    gl.viewport(0, 0, width, height);
  }
}

function convertHexArrToRGBAArr(hexArr) {
  const palette = new Float32Array(hexArr.length * 4);

  for (let i = 0; i < hexArr.length * 4; i += 4) {
    const hexStr = hexArr[i / 4];

    const red = parseInt(hexStr[1], 16) * 16 + parseInt(hexStr[2], 16);
    const green = parseInt(hexStr[3], 16) * 16 + parseInt(hexStr[4], 16);
    const blue = parseInt(hexStr[5], 16) * 16 + parseInt(hexStr[6], 16);

    palette[i] = red / 255;
    palette[i + 1] = green / 255;
    palette[i + 2] = blue / 255;
    palette[i + 3] = 1.0;
  }

  return palette;
}

function maxDistanceChannel(pallete) {
  let result = new Float32Array(3);
  let rChannel = [];
  let gChannel = [];
  let bChannel = [];

  for (let i = 0; i < pallete.length; i += 4) {
    rChannel.push(pallete[i]);
    gChannel.push(pallete[i + 1]);
    bChannel.push(pallete[i + 2]);
  }

  result[0] = maxDistanceInArray(rChannel);
  result[1] = maxDistanceInArray(gChannel);
  result[2] = maxDistanceInArray(bChannel);

  return result;
}

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

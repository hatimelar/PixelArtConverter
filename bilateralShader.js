"use strict";
var vertexShaderSource = `#version 300 es

// an attribute is an input (in) to a vertex shader.
// It will receive data from a buffer
in vec2 a_position;
in vec2 a_texCoord;

// Used to pass in the resolution of the canvas
uniform vec2 u_resolution;
uniform float u_flipY;

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

  gl_Position = vec4(clipSpace * vec2(1, u_flipY), 0, 1);

  // pass the texCoord to the fragment shader
  // The GPU will interpolate this value between points.
  v_texCoord = a_texCoord;
}
`;

var fragmentShaderSource = `#version 300 es
precision highp float;

uniform sampler2D u_image;
uniform float u_kernel_0[25];
uniform float u_kernel_1[256];
uniform int u_vRun;

in vec2 v_texCoord;
out vec4 outColor;

void main() {
  if (u_vRun==-1){
    outColor = texture(u_image, v_texCoord);
  }
  else {
    vec2 onePixel = vec2(1.0) / vec2(textureSize(u_image, 0));
    vec4 pixelColor  =  texture(u_image, v_texCoord);
    vec4 colorSum = vec4(0.0);
    float wR = 0.0;
    float wG = 0.0;
    float wB = 0.0;
    int fail=0;
    for (int i = -16; i <= 16; i++) {
        vec2 neighborCoord;
        if (u_vRun == 1){
          neighborCoord = v_texCoord + onePixel * vec2(0.0, float(i));
          if (neighborCoord.y<0. || neighborCoord.y>1.){
            fail=1;
            break;
          }
        }
        else{
          neighborCoord = v_texCoord + onePixel * vec2(float(i), 0.0);
          if (neighborCoord.x<0. || neighborCoord.x>1.){
            fail=1;
            break;
          }
        } 

        vec4 neighborColor = texture(u_image, neighborCoord);

        float g1 = u_kernel_0[abs(i)];
        
        float diffR = abs(neighborColor.r - pixelColor.r);
        float diffG = abs(neighborColor.g - pixelColor.g);
        float diffB = abs(neighborColor.b - pixelColor.b);

        float cWR = u_kernel_1[int(round(diffR * 255.0))] * g1;
        float cWG = u_kernel_1[int(round(diffG * 255.0))] * g1;
        float cWB = u_kernel_1[int(round(diffB * 255.0))] * g1;

        wR += cWR;
        wG += cWG;
        wB += cWB;

        colorSum += neighborColor * vec4(cWR, cWG, cWB, 1.0);
    }

    if ( fail==0) {
        colorSum /= vec4(wR, wG, wB, 1.0);
        outColor = clamp(colorSum, 0.0, 1.0);
    } else {
        outColor = pixelColor;
    }
  }
}
`;

class BilateralShadingData {
  constructor(
    program,
    positionAttributeLocation,
    texCoordAttributeLocation,
    resolutionLocation,
    imageLocation,
    kernel0Location,
    kernel1Location,
    vRunLocation,
    flipYLocation,
    vao,
    positionBuffer,
    originalImageTexture,
    textures,
    frameBuffers
  ) {
    this.program = program;
    this.positionAttributeLocation = positionAttributeLocation;
    this.texCoordAttributeLocation = texCoordAttributeLocation;
    this.resolutionLocation = resolutionLocation;
    this.imageLocation = imageLocation;
    this.kernel0Location = kernel0Location;
    this.kernel1Location = kernel1Location;
    this.vRunLocation = vRunLocation;
    this.flipYLocation = flipYLocation;
    this.vao = vao;
    this.positionBuffer = positionBuffer;
    this.originalImageTexture = originalImageTexture;
    this.textures = textures;
    this.frameBuffers = frameBuffers;
  }
}

function setupBilateralFilterShadingProgram(gl) {
  // setup GLSL program
  var program = createProgramFromSources(gl, [
    vertexShaderSource,
    fragmentShaderSource,
  ]);

  // look up where the vertex data needs to go.
  var positionAttributeLocation = gl.getAttribLocation(program, "a_position");
  var texCoordAttributeLocation = gl.getAttribLocation(program, "a_texCoord");

  // lookup uniforms
  var resolutionLocation = gl.getUniformLocation(program, "u_resolution");
  var imageLocation = gl.getUniformLocation(program, "u_image");
  var kernel0Location = gl.getUniformLocation(program, "u_kernel_0[0]");
  var kernel1Location = gl.getUniformLocation(program, "u_kernel_1[0]");
  var vRunLocation = gl.getUniformLocation(program, "u_vRun");
  var flipYLocation = gl.getUniformLocation(program, "u_flipY");

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

  // create 2 textures and attach them to framebuffers.
  var textures = [];
  var frameBuffers = [];
  for (var ii = 0; ii < 2; ++ii) {
    var texture = createAndSetupTexture(gl);
    textures.push(texture);

    // Create a framebuffer
    var fbo = gl.createFramebuffer();
    frameBuffers.push(fbo);
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);

    // Attach a texture to it.
    var attachmentPoint = gl.COLOR_ATTACHMENT0;
    gl.framebufferTexture2D(
      gl.FRAMEBUFFER,
      attachmentPoint,
      gl.TEXTURE_2D,
      texture,
      0
    );
  }

  return new BilateralShadingData(
    program,
    positionAttributeLocation,
    texCoordAttributeLocation,
    resolutionLocation,
    imageLocation,
    kernel0Location,
    kernel1Location,
    vRunLocation,
    flipYLocation,
    vao,
    positionBuffer,
    originalImageTexture,
    textures,
    frameBuffers
  );
}

function bilateralSmoothing(
  gl,
  bilateralShadingData,
  image,
  buffer,
  texture,
  iterations
) {
  gl.bindTexture(gl.TEXTURE_2D, bilateralShadingData.originalImageTexture);
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

  for (let i = 0; i < 2; i++) {
    gl.bindTexture(gl.TEXTURE_2D, bilateralShadingData.textures[i]);
    // make the texture the same size as the image
    var mipLevel = 0; // the largest mip
    var internalFormat = gl.RGBA; // format we want in the texture
    var border = 0; // must be 0
    var srcFormat = gl.RGBA; // format of data we are supplying
    var srcType = gl.UNSIGNED_BYTE; // type of data we are supplying
    var data = null; // no data = create a blank texture
    gl.texImage2D(
      gl.TEXTURE_2D,
      mipLevel,
      internalFormat,
      gl.canvas.width,
      gl.canvas.height,
      border,
      srcFormat,
      srcType,
      data
    );
  }
  gl.bindTexture(gl.TEXTURE_2D, texture);
  // make the texture the same size as the image
  var mipLevel = 0; // the largest mip
  var internalFormat = gl.RGBA; // format we want in the texture
  var border = 0; // must be 0
  var srcFormat = gl.RGBA; // format of data we are supplying
  var srcType = gl.UNSIGNED_BYTE; // type of data we are supplying
  var data = null; // no data = create a blank texture
  gl.texImage2D(
    gl.TEXTURE_2D,
    mipLevel,
    internalFormat,
    gl.canvas.width,
    gl.canvas.height,
    border,
    srcFormat,
    srcType,
    data
  );

  // Bind the position buffer so gl.bufferData that will be called
  // in setRectangle puts data in the position buffer
  gl.bindBuffer(gl.ARRAY_BUFFER, bilateralShadingData.positionBuffer);

  // Set a rectangle the same size as the image.
  setRectangle(gl, 0, 0, gl.canvas.width, gl.canvas.height);

  const sigS = 8;
  const sigR = 8;

  const diameter = Math.ceil(sigS * 2) * 2 + 1;
  const center = (diameter - 1) / 2;

  const kernel0Lookup = new Float32Array(center + 1);
  for (let i = 0; i < center + 1; i++) {
    kernel0Lookup[i] = gaussValAt2(i, sigS);
  }
  const kernel1Lookup = new Float32Array(256);
  for (let i = 0; i < 256; i++) kernel1Lookup[i] = gaussValAt2(i, sigR);

  drawEffects();

  function drawEffects() {
    // Tell WebGL how to convert from clip space to pixels
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

    // Clear the canvas
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    // Tell it to use our program (pair of shaders)
    gl.useProgram(bilateralShadingData.program);

    // Bind the attribute/buffer set we want.
    gl.bindVertexArray(bilateralShadingData.vao);

    // start with the original image on unit 0
    gl.activeTexture(gl.TEXTURE0 + 0);
    gl.bindTexture(gl.TEXTURE_2D, bilateralShadingData.originalImageTexture);

    // Tell the shader to get the texture from texture unit 0
    gl.uniform1i(bilateralShadingData.imageLocation, 0);

    // don't y flip images while drawing to the textures
    gl.uniform1f(bilateralShadingData.flipYLocation, 1);

    gl.uniform1fv(bilateralShadingData.kernel0Location, kernel0Lookup);
    gl.uniform1fv(bilateralShadingData.kernel1Location, kernel1Lookup);

    // loop through each effect we want to apply.
    var count = 0;

    for (let i = 0; i < iterations * 2; i++) {
      gl.uniform1i(bilateralShadingData.vRunLocation, count % 2);
      setFramebuffer(
        bilateralShadingData.frameBuffers[count % 2],
        gl.canvas.width,
        gl.canvas.height
      );
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      draw();
      gl.bindTexture(gl.TEXTURE_2D, bilateralShadingData.textures[count % 2]);
      count++;
    }

    setFramebuffer(buffer, gl.canvas.width, gl.canvas.height);

    // Clear the buffer
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.uniform1f(bilateralShadingData.flipYLocation, -1);
    gl.uniform1i(bilateralShadingData.vRunLocation, -1);

    draw();

    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }

  function setFramebuffer(fbo, width, height) {
    // make this the framebuffer we are rendering to.
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);

    // Tell the shader the resolution of the framebuffer.
    gl.uniform2f(bilateralShadingData.resolutionLocation, width, height);

    // Tell WebGL how to convert from clip space to pixels
    gl.viewport(0, 0, width, height);
  }

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
}

function gaussValAt2(x, sig) {
  const exponentX = -(x ** 2) / (2 * sig ** 2);

  return (1 / (sig ** 2 * 2 * Math.PI)) * Math.exp(exponentX);
}

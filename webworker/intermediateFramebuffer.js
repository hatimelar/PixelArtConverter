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

function getIntermediateFramebuffer(gl) {
  var texture = createAndSetupTexture(gl);

  // Create a framebuffer
  var fbo = gl.createFramebuffer();
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
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  return [fbo, texture];
}

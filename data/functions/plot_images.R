# A collection of functions to plot the output SimSpin images with appropriate
# colour bars
# Kate Harborne 09/03/22

library(magicaxis)

.image_nan <- function(z, zlim, col, na.color='gray', ...){
  #' Function for plotting an image with NAs as a set colour
  #' @param z array corresponding to the image to be plotted
  #' @param zlim colour range associated with the image for 
  #              the colour palette
  #' @param col the colour palette to be used for the image
  #' @param na.color string describing the color used for NA
  
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.na <- zlim[2] + zstep # new z for NA
  
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[2] <- zlim[2] + zstep # extend top limit to include the new value na
  
  col <- c(col, na.color) # we construct the new color range by including: na.color
  
  magicaxis::magimage(x=z, zlim=zlim, col=col, ...) # we finally call image(...)
}

magcolbar <-
  function(position='topright',range=c(0,1),orient='v',log=FALSE,col=hcl.colors(21),scale=c(1/4,1/20), inset=1/40,labN=5,title='',titleshift=0,centrealign='rb',clip='',cex=1,...){
    usercoord=par()$usr
    xlogcheck=FALSE;ylogcheck=FALSE
    if(par()$xlog){par(xlog=FALSE);par(usr=c(log10(par()$usr[1:2]),par()$usr[3:4]));xlogcheck=TRUE}
    if(par()$ylog){par(ylog=FALSE);par(usr=c(par()$usr[1:2],log10(par()$usr[3:4])));ylogcheck=TRUE}
    
    xlo=usercoord[1];xhi=usercoord[2];ylo=usercoord[3];yhi=usercoord[4]
    xdiff=xhi-xlo;ydiff=yhi-ylo
    
    if(orient=='h'){
      xl=xlo+xdiff/2-xdiff*scale[1]/2
      yb=ylo+ydiff/2-ydiff*scale[2]/2
      xr=xlo+xdiff/2+xdiff*scale[1]/2
      yt=ylo+ydiff/2+ydiff*scale[2]/2
      align=centrealign
      if(length(grep('bottom',position))>0){align='lt'; yb=ylo+ydiff*inset;yt=ylo+ydiff*inset+ydiff*scale[2]}
      if(length(grep('top',position))>0){align='rb';yb=yhi-ydiff*inset-ydiff*scale[2];yt=yhi-ydiff*inset}
      if(length(grep('left',position))>0){xl=xlo+xdiff*inset;xr=xlo+xdiff*inset+xdiff*scale[1]}
      if(length(grep('right',position))>0){xl=xhi-xdiff*inset-xdiff*scale[1];xr=xhi-xdiff*inset}
    }
    
    if(orient=='v'){
      xl=xlo+xdiff/2-xdiff*scale[2]/2
      yb=ylo+ydiff/2-ydiff*scale[1]/2
      xr=xlo+xdiff/2+xdiff*scale[2]/2
      yt=ylo+ydiff/2+ydiff*scale[1]/2
      align=centrealign
      if(length(grep('bottom',position))>0){yb=ylo+ydiff*inset;yt=ylo+ydiff*inset+ydiff*scale[1]}
      if(length(grep('top',position))>0){yb=yhi-ydiff*inset-ydiff*scale[1];yt=yhi-ydiff*inset}
      if(length(grep('left',position))>0){align='lt';xl=xlo+xdiff*inset;xr=xlo+xdiff*inset+xdiff*scale[2]}
      if(length(grep('right',position))>0){align='rb';xl=xhi-xdiff*inset-xdiff*scale[2];xr=xhi-xdiff*inset}
      #if(length(grep('left',position))>0){align='rb';xl=xlo+xdiff*inset;xr=xlo+xdiff*inset+xdiff*scale[2]}
      #if(length(grep('right',position))>0){align='lt';xl=xhi-xdiff*inset-xdiff*scale[2];xr=xhi-xdiff*inset}
    }
    
    rangetemp=range
    legend=maglab(rangetemp,labN,log=log,trim=F)
    if(log){
      rangetemp=log10(range)
      legend$labat=log10(legend$labat)
    }
    
    roughNscale=(max(legend$labat)-min(legend$labat))/(rangetemp[2]-rangetemp[1])
    if(clip=='bg'){clip='NA'}
    
    colremap=magmap(data=seq(min(legend$labat),max(legend$labat),length=length(col)*roughNscale),locut=rangetemp[1],hicut=rangetemp[2],type='num',range=c(1,length(col)),clip=clip)$map
    col=col[round(colremap,digits=0)]
    
    if(orient=='v'){plotrix::color.legend(xl,yb,xr,yt,legend=legend$exp,rect.col=col,cex=cex,align=align,gradient='y',...)}
    if(orient=='h'){plotrix::color.legend(xl,yb,xr,yt,legend=legend$exp,rect.col=col,cex=cex,align=align,gradient='x',...)}
    
    if(orient=='v' & align=='lt'){text(xl-(1+titleshift)*xdiff/20,(yt+yb)/2,labels=title,adj=c(0.5,0.5),srt=90,cex=cex, xpd=NA)}
    if(orient=='v' & align=='rb'){text(xr+(1+titleshift)*xdiff/20,(yt+yb)/2,labels=title,adj=c(0.5,0.5),srt=-90,cex=cex, xpd=NA)}
    if(orient=='h' & align=='lt'){text((xl+xr)/2,yt+(titleshift)*ydiff/20, labels=title,adj=c(0.5,0.5),srt=0,cex=cex, xpd=NA)}
    if(orient=='h' & align=='rb'){text((xl+xr)/2,yb-(titleshift)*ydiff/20,labels=title,adj=c(0.5,0.5),srt=0,cex=cex, xpd=NA)}
    
    par(xlog=xlogcheck)
    par(ylog=ylogcheck)
    par(usr=usercoord)
  }


plot_flux <- function(flux_image, fig = c(0,1,0,1), new=F,
                      units = expression("r-band Flux, CGS"), zlim=NA, 
                      main="", radii = NA, na.color = "white", legend = T, 
                      titleshift = -4, labN=5, cex=1,...){
  #' Function to plot flux/mass maps with cmocean balance scaling
  #' @param flux_image Numeric array containing the flux/mass image
  #' @param fig Numeric array of length 4 describing the boundary of the image
  #' @param new Boolean. Should the image be added to the existing plot? Default
  #'            is FALSE.
  #' @param units String describing the units of the values contained 
  #'              in the image  
  #' @param radii list - containing a,b,ang (if wishing to plot half-mass radii ellipse)
  
  Flux = flux_image
  im_dim = dim(Flux)%/%2
  flux_val = c(min(Flux, na.rm = T), max(Flux, na.rm = T))
  flux_map_cols = cmocean::cmocean("thermal", version = "2.0", start = .1, end=1)(50)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = Flux, zlim = if(is.na(zlim[1])){flux_val}else{zlim}, 
             col = flux_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd = 3)
  }
  if (legend){
    magcolbar(position = "bottom", range = if(is.na(zlim[1])){flux_val}else{zlim}, scale = c(1, 1/20),
              col = flux_map_cols, orient = "h", inset = -1/20, labN=labN, title = units, cex=cex,
              titleshift = titleshift)
  }  
}

plot_velocity <- function(velocity_image, fig = c(0,1,0,1), new=F,
                          units = expression("velocity"[LOS] * ", km s"^{-1}), main="", 
                          radii = NA, na.color = "white", zlim = NA, legend = T, 
                          titleshift = -4, labN=5, cex=1, ...){
  #' Function to plot velocity maps with cmocean balance scaling
  #' @param velocity_image Numeric matrix containing the velocity image
  #' @param fig Numeric array of length 4 describing the boundary of the image
  #' @param new Boolean. Should the image be added to the existing plot? Default
  #'            is FALSE.
  #' @param units String describing the units of the values contained 
  #'              in the image  
  #' @param radii list - containing a,b,ang (if wishing to plot half-mass radii ellipse)
  
  V = velocity_image  
  im_dim = dim(V)%/%2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))
  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd = 3)
  }
  if (legend){
    magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
              col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units, cex=cex,
              titleshift = titleshift)
  }
  
}

plot_dispersion <- function(dispersion_image, fig = c(0,1,0,1), new=F,
                            units = expression("dispersion"[LOS] * ", km s"^{-1}), main="", 
                            radii = NA, na.color = "white", zlim = NA, legend=T, 
                            titleshift = -4, labN=5, cex=1, ...){
  #' Function to plot dispersion maps with cmocean solar scaling
  #' @param dispersion_image Numeric matrix containing the dispersion image
  #' @param fig Numeric array of length 4 describing the boundary of the image
  #' @param new Boolean. Should the image be added to the existing plot? Default
  #'            is FALSE.
  #' @param units String describing the units of the values contained 
  #'              in the image  
  #' @param radii list - containing a,b,ang (if wishing to plot half-mass radii ellipse)
  #' 
  disp_map = dispersion_image
  im_dim = dim(disp_map)%/%2
  disp_val = c(floor(min(disp_map, na.rm=T)), max(disp_map, na.rm=T))
  disp_map_cols = cmocean::cmocean("solar", version = "2.0", start = .1, end=.9)(100)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = disp_map, zlim = if(is.na(zlim[1])){disp_val}else{zlim}, 
             col = disp_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd = 3)
  }
  if (legend){
    magcolbar(position = "bottom", range = if(is.na(zlim[1])){disp_val}else{zlim}, scale = c(1, 1/20),
              col = disp_map_cols, orient = "h", inset = -1/20, labN=labN, title = units, 
              titleshift = titleshift, cex=cex)
  }
  
}

plot_age <- function(age_image, fig = c(0,1,0,1), new=F,
                     units = expression("Age, Gyr"), main="", radii = NA,
                     na.color = "white", ...){
  #' Function to plot dispersion maps with cmocean solar scaling
  #' @param age_image Numeric matrix containing the dispersion image
  #' @param fig Numeric array of length 4 describing the boundary of the image
  #' @param new Boolean. Should the image be added to the existing plot? Default
  #'            is FALSE.
  #' @param units String describing the units of the values contained 
  #'              in the image  
  
  age_map = age_image
  im_dim = dim(age_map)%/%2
  
  age_val = c(floor(min(age_map, na.rm=T)/2), max(age_map, na.rm=T))
  age_map_cols = cmocean::cmocean("deep", version = "2.0", start = .1, end=1)(50)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = age_map, zlim = age_val, col = age_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd = 3)
  }
  fields::image.plot(legend.only = TRUE, zlim = c(age_val), col = age_map_cols,
                     horizontal = TRUE, family="serif", font=1,
                     legend.lab = units)
  
}

plot_metallicity <- function(metallicity_image, fig = c(0,1,0,1), new=F,
                             units = expression("log10(Z/Z"[solar]*")"), main="", 
                             na.color = "white", ...){
  #' Function to plot dispersion maps with cmocean solar scaling
  #' @param metallicity_image Numeric matrix containing the dispersion image
  #' @param fig Numeric array of length 4 describing the boundary of the image
  #' @param new Boolean. Should the image be added to the existing plot? Default
  #'            is FALSE.
  #' @param units String describing the units of the values contained 
  #'              in the image  
  
  met_map = metallicity_image
  met_val = c(floor(min(met_map, na.rm=T)/2), max(met_map, na.rm=T))
  met_map_cols = cmocean::cmocean("dense", version = "2.0", start = .1, end=1)(50)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = met_map, zlim = met_val, col = met_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  fields::image.plot(legend.only = TRUE, zlim = c(met_val), col = met_map_cols,
                     horizontal = TRUE, family="serif", font=1,
                     legend.lab = units)
  
}

plot_h3   <- function(h3_image, fig = c(0,1,0,1), new=F,
                      units = expression("h"[3]), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, cex=1, ...){
  
  V = h3_image
  im_dim = dim(V)%/%2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))
  
  if (all(vel_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }
  
  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd=3)
  }
  if (legend){
    magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
               col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units, cex=cex,
               titleshift = titleshift)
  }
  
}

plot_h4   <- function(h4_image, fig = c(0,1,0,1), new=F,
                      units = expression("h"[4]), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, cex=1,...){
  
  V = h4_image
  im_dim = dim(V)%/%2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))
  
  if (all(vel_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }
  
  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)
  
  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = "red", density = NULL, lwd = 3)
  }
  if (legend){
    magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
               col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units, cex=cex,
               titleshift = titleshift)
  }
  
}

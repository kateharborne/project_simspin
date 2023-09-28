# Generating the plots for Case Study 1
# Kate Harborne 
#
# Running this file will generate the plots using data for constructed in the  
# Case Study 1 file as shown in the paper (section 3.1 test 1)
#
################################################################################
#
# Installation of packages below may be necessary if you do not have them 
# already installed. Installing reticulate may require further support and an 
# installation of python3 and the numpy package. This can be done using conda. 
# Please refer to the following address for futher help with this:
# https://rstudio.github.io/reticulate/
# Uncomment the lines below to install the necessary modules.

# install.packages("devtools")
# devtools::install_github("kateharborne/SimSpin@v2.6.0")
# install.packages("magicaxis")
# install.packages("reticulate")

################################################################################

# Loading packges -------------------------------------------------------------- 

library(SimSpin)
library(magicaxis)
library(reticulate)
np = import("numpy", convert = T)  

# Sourcing plotting and result functions ---------------------------------------

source("data/functions/plot_images.R")
source("data/functions/plot_array.R")
source("data/functions/summarise_results.R")

# Reference directories --------------------------------------------------------

dir_plots = "data/plots/"
ppxf_out = "data/cubes/ppxf_output/fwhm0_lowz/"
simspin_out = "data/cubes/fwhm0_lowz/"

w_array = plot_array(plot_size = 0.31, n_plots = 4)

# Used obs_analsysis to get an isophote for Re ellipses ------------------------
lowz_scale_EMILES = 0.15; s=2
lowz_scale_BC03 = 0.2125
scale = lowz_scale_EMILES/lowz_scale_BC03

disk_ellipse_data_EMILES = list("a" = 13.9, "b" = 9.7, "ang" = 0)
bulge_ellipse_data_EMILES = list("a" = 12.2, "b" = 10.8, "ang" = 0)
disk_ellipse_data_BC03 = list("a" = 13.9*scale, "b" = 9.7*scale, "ang" = 0)
bulge_ellipse_data_BC03 = list("a" = 12.2*scale, "b" = 10.8*scale, "ang" = 0)

# Case Study 1: (1a) EMILES templates, FWHM = 0, low redshift (DISK) -----------
# disk_lowz_fwhm0_EMILES = summarise_results(fname = "disk_EMILES_spectral_fwhm0_lowz", 
#                                            folder = "disk_EMILES/")
# 

disk_lowz_fwhm0_EMILES = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm0_EMILES.Rdata"))

# Case Study 1: (1b) BC03hr templates, FWHM = 0, low redshift (DISK) -----------
# disk_lowz_fwhm0_BC03hr = summarise_results(fname = "disk_BC03hr_spectral_fwhm0_lowz", 
#                                            folder = "disk_BC03hr/")
# 
# saveRDS(disk_lowz_fwhm0_BC03hr, paste0(simspin_out, "summary_disk_lowz_fwhm0_BC03hr.Rdata"))

disk_lowz_fwhm0_BC03hr = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm0_BC03hr.Rdata"))

# Case Study 1: (2a) EMILES templates, FWHM = 0, low redshift (BULGE) ----------
# bulge_lowz_fwhm0_EMILES = summarise_results(fname = "bulge_EMILES_spectral_fwhm0_lowz", 
#                                             folder = "bulge_EMILES/")
# 
# saveRDS(bulge_lowz_fwhm0_EMILES, paste0(simspin_out, "summary_bulge_lowz_fwhm0_EMILES.Rdata"))

bulge_lowz_fwhm0_EMILES = readRDS(paste0(simspin_out, "summary_bulge_lowz_fwhm0_EMILES.Rdata"))

# Case Study 1: (2b) BC03hr templates, FWHM = 0, low redshift (BULGE) ----------
# bulge_lowz_fwhm0_BC03hr = summarise_results(fname = "bulge_BC03hr_spectral_fwhm0_lowz", 
#                                             folder = "bulge_BC03hr/")
# 
# saveRDS(bulge_lowz_fwhm0_BC03hr, paste0(simspin_out, "summary_bulge_lowz_fwhm0_BC03hr.Rdata"))

bulge_lowz_fwhm0_BC03hr = readRDS(paste0(simspin_out, "summary_bulge_lowz_fwhm0_BC03hr.Rdata"))

# Case Study 1: (3a) EMILES templates, FWHM = 0, low redshift (BULGE) ----------
# old_bulge_lowz_fwhm0_EMILES = summarise_results(fname = "old_bulge_EMILES_spectral_fwhm0_lowz",
#                                                 folder = "old_bulge_EMILES/")
# 
# saveRDS(old_bulge_lowz_fwhm0_EMILES, paste0(simspin_out, "summary_old_bulge_lowz_fwhm0_EMILES.Rdata"))

old_bulge_lowz_fwhm0_EMILES = readRDS(paste0(simspin_out, "summary_old_bulge_lowz_fwhm0_EMILES.Rdata"))

# Case Study 1: (3b) BC03hr templates, FWHM = 0, low redshift (BULGE) ----------
# old_bulge_lowz_fwhm0_BC03hr = summarise_results(fname = "old_bulge_BC03hr_spectral_fwhm0_lowz",
#                                                 folder = "old_bulge_BC03hr/")
# 
# saveRDS(old_bulge_lowz_fwhm0_BC03hr, paste0(simspin_out, "summary_old_bulge_lowz_fwhm0_BC03hr.Rdata"))

old_bulge_lowz_fwhm0_BC03hr = readRDS(paste0(simspin_out, "summary_old_bulge_lowz_fwhm0_BC03hr.Rdata"))

# Maps =========================================================================
# disk EMILES ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_disk_velocities_lowz_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))
# velocity
plot_velocity(disk_lowz_fwhm0_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = disk_lowz_fwhm0_EMILES$vel_lims, labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_EMILES, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_EMILES))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm0_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = disk_lowz_fwhm0_EMILES$vel_lims, labN = 2, titleshift = -5, units="",#
              radii = disk_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm0_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm0_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-10,10),
              radii = disk_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -65, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 305, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm0_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm0_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES, cex=1.2)

plot_velocity(disk_lowz_fwhm0_EMILES$res_dispersion, zlim = c(-20,20),
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -90, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm0_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#disk_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm0_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#zlim = disk_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES, cex=1.2)

plot_h3(disk_lowz_fwhm0_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, units="",#
        #units = expression("residual, h"[3]),
        zlim = c(-0.1,0.1),#
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_EMILES$res_h3, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.1,0.1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -75, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm0_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm0_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm0_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm0_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES, cex=1.2)

plot_h4(disk_lowz_fwhm0_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        zlim = c(-0.1,0.1),
        #units = expression("residual, h"[4]),
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_EMILES$res_h4, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.1,0.1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -65, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()

# disk BC03hr ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_disk_velocities_lowz_BC03.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

# Velocity
plot_velocity(disk_lowz_fwhm0_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_BC03, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_BC03)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_BC03))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm0_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units = "",
              radii = disk_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm0_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm0_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, #units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, units = "",
              zlim = c(-50,50),
              radii = disk_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-2, 2), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 145, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# Dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm0_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm0_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03, cex=1.2)

plot_velocity(disk_lowz_fwhm0_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, 
              units="",#units = expression("residual, km s"^{-1}),
              zlim = c(-50, 50),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -18, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm0_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.06,0.06), labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm0_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = disk_lowz_fwhm0_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03, cex=1.2)

plot_h3(disk_lowz_fwhm0_BC03hr$res_h3,
        zlim = disk_lowz_fwhm0_BC03hr$h3_lims,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",#units = expression("residual, h"[3]), 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_BC03hr$res_h3,  yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -22, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

#h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm0_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.15,0.15), labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm0_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = disk_lowz_fwhm0_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03, cex=1.2)

plot_h4(disk_lowz_fwhm0_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",#units = expression("residual, h"[4]),
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm0_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -19, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()

# bulge EMILES -----------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_bulge_velocities_lowz_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

plot_velocity(bulge_lowz_fwhm0_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-40,40), units="",#bulge_lowz_fwhm0_EMILES$vel_lims, 
              labN = 2, titleshift = -5,
              radii = bulge_ellipse_data_EMILES, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_EMILES))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(bulge_lowz_fwhm0_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-40,40), units="",
              labN = 2, titleshift = -5,
              radii = bulge_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(bulge_lowz_fwhm0_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(bulge_lowz_fwhm0_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units = "",#expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-40,40),
              radii = bulge_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1,1), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
text(0, 330, labels="Residual",adj=c(0.5,0.5), xpd=NA)


par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(bulge_lowz_fwhm0_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = bulge_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units="",
                radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(bulge_lowz_fwhm0_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = bulge_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units="",
                radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_velocity(bulge_lowz_fwhm0_EMILES$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",
              labN = 2, titleshift = -5,
              radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)


par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(bulge_lowz_fwhm0_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#bulge_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5,  units="",
        radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(bulge_lowz_fwhm0_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#bulge_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5,  units="",
        radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_h3(bulge_lowz_fwhm0_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units="",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5,
        radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_EMILES$res_h3, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)

par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(bulge_lowz_fwhm0_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = bulge_lowz_fwhm0_EMILES$h4_lims, labN = 2, titleshift = -5, units="",
        radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(bulge_lowz_fwhm0_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = bulge_lowz_fwhm0_EMILES$h4_lims, labN = 2, titleshift = -5, units="",
        radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_h4(bulge_lowz_fwhm0_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5,
        radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_EMILES$res_h4, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)

dev.off()


# bulge BC03hr -----------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_bulge_velocities_lowz_BC03hr.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

plot_velocity(bulge_lowz_fwhm0_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-40,40),#bulge_lowz_fwhm0_BC03hr$vel_lims, 
              labN = 2, titleshift = -5, units="",
              radii=bulge_ellipse_data_BC03, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_BC03)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_BC03))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(bulge_lowz_fwhm0_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-40,40),
              labN = 2, titleshift = -5, units="",
              radii=bulge_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(bulge_lowz_fwhm0_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(bulge_lowz_fwhm0_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",
              labN = 2, titleshift = -5, zlim = c(-40,40),
              radii=bulge_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 80, labels="Residual",adj=c(0.5,0.5), xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)

par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(bulge_lowz_fwhm0_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = bulge_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(bulge_lowz_fwhm0_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = bulge_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii=bulge_ellipse_data_BC03, cex=1.2)

plot_velocity(bulge_lowz_fwhm0_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",
              labN = 2, titleshift = -5, zlim = c(-50,50),
              radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size,
        col = "#785EF0", border = NA, xlim = c(-1, 1), yaxt="n",
        mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -15, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)



par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(bulge_lowz_fwhm0_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#bulge_lowz_fwhm0_BC03hr$h3_lims, 
        labN = 2, titleshift = -5, units="",
        radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(bulge_lowz_fwhm0_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#bulge_lowz_fwhm0_BC03hr$h3_lims, 
        labN = 2, titleshift = -5, units="",
        radii=bulge_ellipse_data_BC03, cex=1.2)

plot_h3(bulge_lowz_fwhm0_BC03hr$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units="",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5,
        radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_BC03hr$res_h3, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -15, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)


par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(bulge_lowz_fwhm0_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = bulge_lowz_fwhm0_BC03hr$h4_lims, labN = 2, titleshift = -5, units="",
        radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(bulge_lowz_fwhm0_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = bulge_lowz_fwhm0_BC03hr$h4_lims, labN = 2, titleshift = -5, units="",
        radii=bulge_ellipse_data_BC03, cex=1.2)

plot_h4(bulge_lowz_fwhm0_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",
        zlim = bulge_lowz_fwhm0_BC03hr$h4_lims,
        labN = 2, titleshift = -5,
        radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(bulge_lowz_fwhm0_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -15, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()


# old bulge EMILES -----------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_old_bulge_velocities_lowz_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

plot_velocity(old_bulge_lowz_fwhm0_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-40,40),#old_bulge_lowz_fwhm0_EMILES$vel_lims, 
              labN = 2, titleshift = -5, units = "",
              radii = bulge_ellipse_data_EMILES, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_EMILES))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(old_bulge_lowz_fwhm0_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-40,40), units = "",
              labN = 2, titleshift = -5, 
              radii = bulge_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)

plot_velocity(old_bulge_lowz_fwhm0_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units = "",
              labN = 2, titleshift = -5, zlim = c(-40,40), 
              radii = bulge_ellipse_data_EMILES, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)
text(0, 300, labels="Residual",adj=c(0.5,0.5), xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)



par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(old_bulge_lowz_fwhm0_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = old_bulge_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units = "",
                radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(old_bulge_lowz_fwhm0_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = old_bulge_lowz_fwhm0_EMILES$disp_lims, labN = 2, titleshift = -5, units = "",
                radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_dispersion(old_bulge_lowz_fwhm0_EMILES$res_dispersion,
                fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units = "",
                labN = 2, titleshift = -5, 
                radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)



par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(old_bulge_lowz_fwhm0_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#old_bulge_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5, units = "",
        radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(old_bulge_lowz_fwhm0_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#old_bulge_lowz_fwhm0_EMILES$h3_lims, 
        labN = 2, titleshift = -5, units = "",
        radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_h3(old_bulge_lowz_fwhm0_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5, 
        radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_EMILES$res_h3, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)


par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(old_bulge_lowz_fwhm0_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = old_bulge_lowz_fwhm0_EMILES$h4_lims, labN = 2, titleshift = -5, units = "",
        radii = bulge_ellipse_data_EMILES, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(old_bulge_lowz_fwhm0_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = old_bulge_lowz_fwhm0_EMILES$h4_lims, labN = 2, titleshift = -5, units = "",
        radii = bulge_ellipse_data_EMILES, cex=1.2)

plot_h4(old_bulge_lowz_fwhm0_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units = "",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5, 
        radii = bulge_ellipse_data_EMILES, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_EMILES$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)

dev.off()



# old bulge BC03hr -----------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_1/cs1_old_bulge_velocities_lowz_BC03hr.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

plot_velocity(old_bulge_lowz_fwhm0_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-40,40),#old_bulge_lowz_fwhm0_BC03hr$vel_lims, 
              labN = 2, titleshift = -5, units = "",
              radii=bulge_ellipse_data_BC03, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_BC03)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_BC03))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(old_bulge_lowz_fwhm0_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-40,40),
              labN = 2, titleshift = -5, units = "",
              radii=bulge_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(old_bulge_lowz_fwhm0_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(old_bulge_lowz_fwhm0_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units = "",
              labN = 2, titleshift = -5, zlim = c(-40,40),
              radii=bulge_ellipse_data_BC03, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 150, labels="Residual",adj=c(0.5,0.5), xpd=NA)




par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(old_bulge_lowz_fwhm0_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = old_bulge_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units = "",
                radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(old_bulge_lowz_fwhm0_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = old_bulge_lowz_fwhm0_BC03hr$disp_lims, labN = 2, titleshift = -5, units = "",
                radii=bulge_ellipse_data_BC03, cex=1.2)

plot_velocity(old_bulge_lowz_fwhm0_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units = "",
              labN = 2, titleshift = -5,
              radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), 
        mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)



par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(old_bulge_lowz_fwhm0_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#old_bulge_lowz_fwhm0_BC03hr$h3_lims, 
        labN = 2, titleshift = -5, units = "",
        radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(old_bulge_lowz_fwhm0_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = c(-0.4,0.4),#old_bulge_lowz_fwhm0_BC03hr$h3_lims, 
        labN = 2, titleshift = -5, units = "",
        radii=bulge_ellipse_data_BC03, cex=1.2)

plot_h3(old_bulge_lowz_fwhm0_BC03hr$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",
        zlim = c(-0.4,0.4),
        labN = 2, titleshift = -5,
        radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_BC03hr$res_h3, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)


par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(old_bulge_lowz_fwhm0_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = old_bulge_lowz_fwhm0_BC03hr$h4_lims, labN = 2, titleshift = -5, units = "",
        radii=bulge_ellipse_data_BC03, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(old_bulge_lowz_fwhm0_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = old_bulge_lowz_fwhm0_BC03hr$h4_lims, labN = 2, titleshift = -5, units = "",
        radii=bulge_ellipse_data_BC03, cex=1.2)

plot_h4(old_bulge_lowz_fwhm0_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units = "",
        zlim = old_bulge_lowz_fwhm0_BC03hr$h4_lims,
        labN = 2, titleshift = -5,
        radii=bulge_ellipse_data_BC03, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(old_bulge_lowz_fwhm0_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
abline(v=0)

dev.off()


# histogram summary ------------------------------------------------------------
plot_dist = function(data_for_hist, lty, col){
  h = maghist(data_for_hist, breaks = seq(-1, 1, length.out=50), plot = F)
  lines(h$mids, h$counts, lwd=2, lty=lty, col=col)
}

w_array = plot_array(plot_size = 0.55, n_plots = 2)
h_array = plot_array(plot_size = 0.54, n_plots = 2)

jpeg(filename = paste0(dir_plots, "case_study_1/cs1_histograms.jpeg"), 
     width = 610, height = 650, res = 100)
par(mar = c(4.5,3.5,1.5,0.5))

par(fig = c(w_array[1,], h_array[2,]))
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,220),
        ylab = "Pixels", xlab=expression("residual v"["LOS"]*"/"*paste(Delta)*paste(nu)))
plot_dist(disk_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, 
          col = "#648FFF", lty=1)
plot_dist(disk_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, 
          col = "#648FFF", lty=3)
plot_dist(bulge_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, 
          col = "#DC267F", lty=1)
plot_dist(bulge_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, 
          col = "#DC267F", lty=3)
plot_dist(old_bulge_lowz_fwhm0_EMILES$res_velocity/EMILES_obs$vbin_size, 
          col = "#FFB000", lty=1)
plot_dist(old_bulge_lowz_fwhm0_BC03hr$res_velocity/BC03_obs$vbin_size, 
          col = "#FFB000", lty=3)
legend("topleft", legend = c("Disc", "Bulge", "Old Bulge"),
       col = c("#648FFF", "#DC267F", "#FFB000"), pch=15, bty = "n")
abline(v=0, lwd=2)

par(fig = c(w_array[2,], h_array[2,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,220), yaxt="n",
        #ylab = "Pixels",
        xlab=expression("residual "*paste(sigma)["LOS"]*"/"*paste(Delta)*paste(nu))) 
plot_dist(disk_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, 
          col="#648FFF", lty=1)
plot_dist(disk_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size, 
          col="#648FFF", lty=3)
plot_dist(bulge_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, 
          col = "#DC267F", lty=1)
plot_dist(bulge_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size, 
          col="#DC267F", lty=3)
plot_dist(old_bulge_lowz_fwhm0_EMILES$res_dispersion/EMILES_obs$vbin_size, 
          col = "#FFB000", lty=1)
plot_dist(old_bulge_lowz_fwhm0_BC03hr$res_dispersion/BC03_obs$vbin_size, 
          col="#FFB000", lty=3)
abline(v=0, lwd=2)

legend("topright", legend = c("BC03hr", "EMILES"),
       col = c("black", "black"), lty=c(3,1), lwd=2, bty = "n")

par(fig = c(w_array[1,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,500),
        ylab = "Pixels",
        xlab=expression("residual h"[3]))
plot_dist(disk_lowz_fwhm0_EMILES$res_h3,
          col = "#648FFF", lty=1)
plot_dist(disk_lowz_fwhm0_BC03hr$res_h3,
          col = "#648FFF", lty=3)
plot_dist(bulge_lowz_fwhm0_EMILES$res_h3,
          col = "#DC267F", lty=1)
plot_dist(bulge_lowz_fwhm0_BC03hr$res_h3,
          col = "#DC267F", lty=3)
plot_dist(old_bulge_lowz_fwhm0_EMILES$res_h3,
          col = "#FFB000", lty=1)
plot_dist(old_bulge_lowz_fwhm0_BC03hr$res_h3,
          col = "#FFB000", lty=3)
abline(v=0, lwd=2)

par(fig = c(w_array[2,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,500), yaxt="n",
        #ylab = "Pixels", 
        xlab=expression("residual h"[4]))
plot_dist(disk_lowz_fwhm0_EMILES$res_h4,
          col = "#648FFF", lty=1)
plot_dist(disk_lowz_fwhm0_BC03hr$res_h4,
          col = "#648FFF", lty=3)
plot_dist(bulge_lowz_fwhm0_EMILES$res_h4,
          col = "#DC267F", lty=1)
plot_dist(bulge_lowz_fwhm0_BC03hr$res_h4,
          col = "#DC267F", lty=3)
plot_dist(old_bulge_lowz_fwhm0_EMILES$res_h4,
          col = "#FFB000", lty=1)
plot_dist(old_bulge_lowz_fwhm0_BC03hr$res_h4,
          col = "#FFB000", lty=3)
abline(v=0, lwd=2)

dev.off()
# per particle/chi/res ---------------------------------------------------------
w_array = plot_array(plot_size = 0.41, n_plots = 3)

jpeg(filename = paste0(dir_plots, "case_study_1/cs1_dispersion_residuals.jpeg"), 
     width = 150, height = 250, units = "mm", res=200)

par(fig = c(0,1, w_array[3,]))
magplot(disk_lowz_fwhm0_EMILES$r_bin_labels, disk_lowz_fwhm0_EMILES$res_r_dispersion_med,
        xlab = c("Radius, kpc"), 
        ylab = expression("Residual "*paste(sigma)["LOS"]),
        type = "l", lty=1, lwd = 2, ylim = c(-50, 50),
        col = "#648FFF")
lines(disk_lowz_fwhm0_BC03hr$r_bin_labels, disk_lowz_fwhm0_BC03hr$res_r_dispersion_med,
      lty=2, lwd = 2,
      col = "#648FFF")
lines(bulge_lowz_fwhm0_EMILES$r_bin_labels, bulge_lowz_fwhm0_EMILES$res_r_dispersion_med,
      lty=1, lwd = 2,
      col = "#DC267F")
lines(bulge_lowz_fwhm0_BC03hr$r_bin_labels, bulge_lowz_fwhm0_BC03hr$res_r_dispersion_med,
      lty=2, lwd = 2,
      col = "#DC267F")
lines(old_bulge_lowz_fwhm0_EMILES$r_bin_labels, old_bulge_lowz_fwhm0_EMILES$res_r_dispersion_med,
      lty=1, lwd = 2,
      col = "#FFB000")
lines(old_bulge_lowz_fwhm0_BC03hr$r_bin_labels, old_bulge_lowz_fwhm0_BC03hr$res_r_dispersion_med,
      lty=2, lwd = 2,
      col = "#FFB000")
abline(h=0, lwd=1)
legend("topleft", legend = c("Disc", "Bulge", "Old Bulge"),
       col = c("#648FFF", "#DC267F", "#FFB000"), pch=15, bty = "n")
legend("topright", legend = c("BC03hr", "EMILES"),
       col = c("black", "black"), lty=c(2,1), bty = "n")

par(fig = c(0,1, w_array[2,]), new=T)
magplot(disk_lowz_fwhm0_EMILES$n_bin_labels, disk_lowz_fwhm0_EMILES$res_n_dispersion_med,
        xlab = c("Number of particles per pixel"), 
        ylab = expression("Residual "*paste(sigma)["LOS"]),
        type = "l", lty=1, lwd = 2, ylim = c(-50, 50), xlim = c(375, 2750),
        col = "#648FFF")
lines(disk_lowz_fwhm0_BC03hr$n_bin_labels, disk_lowz_fwhm0_BC03hr$res_n_dispersion_med,
      lty=2, lwd = 2,
      col = "#648FFF")
lines(bulge_lowz_fwhm0_EMILES$n_bin_labels, bulge_lowz_fwhm0_EMILES$res_n_dispersion_med,
      lty=1, lwd = 2,
      col = "#DC267F")
lines(bulge_lowz_fwhm0_BC03hr$n_bin_labels, bulge_lowz_fwhm0_BC03hr$res_n_dispersion_med,
      lty=2, lwd = 2,
      col = "#DC267F")
lines(old_bulge_lowz_fwhm0_EMILES$n_bin_labels, old_bulge_lowz_fwhm0_EMILES$res_n_dispersion_med,
      lty=1, lwd = 2,
      col = "#FFB000")
lines(old_bulge_lowz_fwhm0_BC03hr$n_bin_labels, old_bulge_lowz_fwhm0_BC03hr$res_n_dispersion_med,
      lty=2, lwd = 2,
      col = "#FFB000")
abline(h=0, lwd=1)

par(fig = c(0,1, w_array[1,]), new=T)
magplot(disk_lowz_fwhm0_EMILES$c_bin_labels, disk_lowz_fwhm0_EMILES$res_c_dispersion_med,
        xlab = expression(paste(chi)^2*" / DOF"), 
        ylab = expression("Residual "*paste(sigma)["LOS"]),
        type = "l", lty=1, lwd = 2, ylim = c(-50, 50), xlim = c(1, 10),
        col = "#648FFF")
lines(disk_lowz_fwhm0_BC03hr$c_bin_labels, disk_lowz_fwhm0_BC03hr$res_c_dispersion_med,
      lty=2, lwd = 2,
      col = "#648FFF")
lines(bulge_lowz_fwhm0_EMILES$c_bin_labels, bulge_lowz_fwhm0_EMILES$res_c_dispersion_med,
      lty=1, lwd = 2,
      col = "#DC267F")
lines(bulge_lowz_fwhm0_BC03hr$c_bin_labels, bulge_lowz_fwhm0_BC03hr$res_c_dispersion_med,
      lty=2, lwd = 2,
      col = "#DC267F")
lines(old_bulge_lowz_fwhm0_EMILES$c_bin_labels, old_bulge_lowz_fwhm0_EMILES$res_c_dispersion_med,
      lty=1, lwd = 2,
      col = "#FFB000")
lines(old_bulge_lowz_fwhm0_BC03hr$c_bin_labels, old_bulge_lowz_fwhm0_BC03hr$res_c_dispersion_med,
      lty=2, lwd = 2,
      col = "#FFB000")
abline(h=0, lwd=1)

dev.off()

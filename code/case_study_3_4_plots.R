# Generating the plots for Case Study 3 and 4
# Kate Harborne 
#
# Running this file will generate the plots using data for constructed in the  
# Case Study 3 and 4 file as shown in the paper (section 3.1 test 3 and 4)
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

w_array = plot_array(plot_size = 0.31, n_plots = 4)

# Getting vscale  --------------------------------------------------------------
SAMI = telescope(type = "IFU", 
                 signal_to_noise = 30, 
                 lsf_fwhm = 3.61, 
                 fov = 15, 
                 wave_range = c(3750,5750))

SAMI_bc = telescope(type = "IFU", 
                    signal_to_noise = 30, 
                    lsf_fwhm = 4.56, 
                    spatial_res = 0.7, 
                    wave_res = 3,
                    fov = 17, 
                    wave_range = c(3750,5750),
                    aperture_shape = "hexagonal")

strategy_one  = observing_strategy(dist_kpc_per_arcsec = 0.3, 
                                   inc_deg = 60, 
                                   blur=F)

EMILES_obs_lowz = observation(SAMI, strategy_one, "velocity")
BC03_obs_lowz = observation(SAMI_bc, strategy_one, "velocity")

strategy_two  = observing_strategy(dist_z = 0.3, 
                                   inc_deg = 60, 
                                   blur=F)

EMILES_obs_highz = observation(SAMI, strategy_two, "velocity")
BC03_obs_highz = observation(SAMI_bc, strategy_two, "velocity")

# Used obs_analsysis to get an isophote for Re ellipses ------------------------
lowz_scale_EMILES = 0.15; s=2
lowz_scale_BC03 = 0.2125
scale = lowz_scale_EMILES/lowz_scale_BC03

disk_ellipse_data_EMILES_lowz = list("a" = 13.9, "b" = 9.7, "ang" = 0)
disk_ellipse_data_BC03_lowz = list("a" = 13.9*scale, "b" = 9.7*scale, "ang" = 0)

# Used obs_analsysis to get an isophote for Re ellipses ------------------------
highz_scale_EMILES = 2.278713; s=2
highz_scale_BC03 = 3.228177
scale = highz_scale_EMILES/highz_scale_BC03

disk_ellipse_data_EMILES_highz = list("a" = 2.8, "b" = 1.4, "ang" = 0)
disk_ellipse_data_BC03_highz = list("a" = 2.8*scale, "b" = 1.4*scale, "ang" = 0)

# Case Study 3: (1a) EMILES templates, FWHM = 3.61, low redshift (DISK) --------

ppxf_out = "data/cubes/ppxf_output/fwhm_lowz/"
simspin_out = "data/cubes/fwhm_lowz/"

# disk_lowz_fwhm_EMILES = summarise_results(fname = "disk_EMILES_spectral_fwhm_lowz", 
#                                            folder = "disk_EMILES/")
# 
# saveRDS(disk_lowz_fwhm_EMILES, paste0(simspin_out, "summary_disk_lowz_fwhm_EMILES.Rdata"))

disk_lowz_fwhm_EMILES = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm_EMILES.Rdata"))

# Case Study 3: (1b) BC03hr templates, FWHM = 4.56, low redshift (DISK) -----
# disk_lowz_fwhm_BC03hr = summarise_results(fname = "disk_BC03hr_spectral_fwhm_lowz", 
#                                            folder = "disk_BC03hr/")
# 
# saveRDS(disk_lowz_fwhm_BC03hr, paste0(simspin_out, "summary_disk_lowz_fwhm_BC03hr.Rdata"))

disk_lowz_fwhm_BC03hr = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm_BC03hr.Rdata"))

# Case Study 3: (2a) EMILES templates, FWHM = 3.61, high redshift (DISK) -----

ppxf_out = "data/cubes/ppxf_output/fwhm_highz/"
simspin_out = "data/cubes/fwhm_highz/"

# disk_highz_fwhm_EMILES = summarise_results(fname = "disk_EMILES_spectral_fwhm_highz", 
#                                           folder = "disk_EMILES/")
# 
# saveRDS(disk_highz_fwhm_EMILES, paste0(simspin_out, "summary_disk_highz_fwhm_EMILES.Rdata"))

disk_highz_fwhm_EMILES = readRDS(paste0(simspin_out, "summary_disk_highz_fwhm_EMILES.Rdata"))

# Case Study 3: (2b) BC03hr templates, FWHM = 4.56, high redshift (DISK) -----
# disk_highz_fwhm_BC03hr = summarise_results(fname = "disk_BC03hr_spectral_fwhm_highz", 
#                                           folder = "disk_BC03hr/")
# 
# saveRDS(disk_highz_fwhm_BC03hr, paste0(simspin_out, "summary_disk_highz_fwhm_BC03hr.Rdata"))
simspin_out = "project_caro/simspin_files/cubes/jun16_2023/fwhm_highz/"
disk_highz_fwhm_BC03hr = readRDS(paste0(simspin_out, "summary_disk_highz_fwhm_BC03hr.Rdata"))

# Case Study 4: (1a) EMILES templates, FWHM = 3.61, low redshift (DISK) -----
ppxf_out = "data/cubes/ppxf_output/fwhm_lowz_blur/"
simspin_out = "data/cubes/fwhm_lowz_blur/"

# disk_lowz_fwhm_blur_EMILES = summarise_results(fname = "disk_EMILES_spectral_fwhm_lowz_blur", 
#                                           folder = "disk_EMILES/")
# 
# saveRDS(disk_lowz_fwhm_blur_EMILES, paste0(simspin_out, "summary_disk_lowz_fwhm_blur_EMILES.Rdata"))

disk_lowz_fwhm_blur_EMILES = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm_blur_EMILES.Rdata"))

# Case Study 4: (1b) BC03hr templates, FWHM = 4.56, low redshift (DISK) -----
# disk_lowz_fwhm_blur_BC03hr = summarise_results(fname = "disk_BC03hr_spectral_fwhm_lowz_blur", 
#                                           folder = "disk_BC03hr/")
# 
# saveRDS(disk_lowz_fwhm_blur_BC03hr, paste0(simspin_out, "summary_disk_lowz_fwhm_blur_BC03hr.Rdata"))

disk_lowz_fwhm_blur_BC03hr = readRDS(paste0(simspin_out, "summary_disk_lowz_fwhm_blur_BC03hr.Rdata"))

# Case Study 4: (2a) EMILES templates, FWHM = 3.61, high redshift (DISK) -----
ppxf_out = "data/cubes/ppxf_output/fwhm_highz_blur/"
simspin_out = "data/cubes/fwhm_highz_blur/"
# disk_highz_fwhm_blur_EMILES = summarise_results(fname = "disk_EMILES_spectral_fwhm_highz_blur", 
#                                            folder = "disk_EMILES/")
# 
# saveRDS(disk_highz_fwhm_blur_EMILES, paste0(simspin_out, "summary_disk_highz_fwhm_blur_EMILES.Rdata"))

disk_highz_fwhm_blur_EMILES = readRDS(paste0(simspin_out, "summary_disk_highz_fwhm_blur_EMILES.Rdata"))

# Case Study 4: (2b) BC03hr templates, FWHM = 4.56, high redshift (DISK) -----
# disk_highz_fwhm_blur_BC03hr = summarise_results(fname = "disk_BC03hr_spectral_fwhm_highz_blur", 
#                                            folder = "disk_BC03hr/")
# 
# saveRDS(disk_highz_fwhm_blur_BC03hr, paste0(simspin_out, "summary_disk_highz_fwhm_blur_BC03hr.Rdata"))

disk_highz_fwhm_blur_BC03hr = readRDS(paste0(simspin_out, "summary_disk_highz_fwhm_blur_BC03hr.Rdata"))

# Maps =========================================================================
# fwhm lowz disk EMILES ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_3/cs3_disk_velocities_lowz_fwhm_EMILES.jpeg"),  
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))
# velocity
plot_velocity(disk_lowz_fwhm_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = disk_lowz_fwhm_EMILES$vel_lims, labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_EMILES))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = disk_lowz_fwhm_EMILES$vel_lims, labN = 2, titleshift = -5, units="",#
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-10,10),
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_EMILES$res_velocity/EMILES_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -65, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 450, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_velocity(disk_lowz_fwhm_EMILES$res_dispersion, zlim = c(-20,20),
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_EMILES$res_dispersion/EMILES_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -90, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#zlim = disk_lowz_fwhm_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_h3(disk_lowz_fwhm_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, units="",#
        #units = expression("residual, h"[3]),
        zlim = c(-0.1,0.1),#
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_EMILES$res_h3, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.2,0.2), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -75, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_h4(disk_lowz_fwhm_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        zlim = c(-0.1,0.1),
        #units = expression("residual, h"[4]),
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_EMILES$res_h4, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.2,0.2), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -65, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()

# fwhm lowz disk BC03hr ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_3/cs3_disk_velocities_lowz_fwhm_BC03.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

# Velocity
plot_velocity(disk_lowz_fwhm_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_BC03)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_BC03))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units = "",
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, #units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, units = "",
              zlim = c(-50,50),
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_BC03hr$res_velocity/BC03_obs_lowz$vbin_size, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 133, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# Dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_BC03hr$disp_lims, labN = 1, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_BC03hr$disp_lims, labN = 1, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_velocity(disk_lowz_fwhm_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, 
              units="",#units = expression("residual, km s"^{-1}),
              zlim = c(-100, 100),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_BC03hr$res_dispersion/BC03_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -18, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.06,0.06), labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = disk_lowz_fwhm_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_h3(disk_lowz_fwhm_BC03hr$res_h3,
        zlim = disk_lowz_fwhm_BC03hr$h3_lims,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",#units = expression("residual, h"[3]), 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_BC03hr$res_h3,  yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -22, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

#h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.15,0.15), labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = disk_lowz_fwhm_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_h4(disk_lowz_fwhm_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",#units = expression("residual, h"[4]),
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -19, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()


# fwhm highz disk EMILES ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_3/cs3_disk_velocities_highz_fwhm_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))
# velocity
plot_velocity(disk_highz_fwhm_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
lines(x = c(s+2, s+2+(10/highz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+4+(10/highz_scale_EMILES))/2, 5, labels="10 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_highz_fwhm_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",#
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(disk_highz_fwhm_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_highz_fwhm_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-10,10),
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_EMILES$res_velocity/EMILES_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -65, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 340, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_highz_fwhm_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_highz_fwhm_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_velocity(disk_highz_fwhm_EMILES$res_dispersion, zlim = c(-20,20),
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_EMILES$res_dispersion/EMILES_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -90, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_highz_fwhm_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#disk_highz_fwhm_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_highz_fwhm_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#zlim = disk_highz_fwhm_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_h3(disk_highz_fwhm_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, units="",#
        #units = expression("residual, h"[3]),
        zlim = c(-0.1,0.1),#
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_EMILES$res_h3, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -75, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_highz_fwhm_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_highz_fwhm_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_highz_fwhm_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_highz_fwhm_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_h4(disk_highz_fwhm_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        zlim = c(-0.1,0.1),
        #units = expression("residual, h"[4]),
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_EMILES$res_h4, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -65, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()



# fwhm highz disk BC03hr ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_3/cs3_disk_velocities_highz_fwhm_BC03.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

# Velocity
plot_velocity(disk_highz_fwhm_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
lines(x = c(s, s+(10/highz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(10/highz_scale_EMILES))/2 + 1, 5, labels="10 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_highz_fwhm_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units = "",
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(disk_highz_fwhm_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_highz_fwhm_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, #units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, units = "",
              zlim = c(-50,50),
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_BC03hr$res_velocity/BC03_obs_highz$vbin_size, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 120, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# Dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_highz_fwhm_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_highz_fwhm_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_velocity(disk_highz_fwhm_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, 
              units="",#units = expression("residual, km s"^{-1}),
              zlim = c(-100, 100),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_BC03hr$res_dispersion/BC03_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -18, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_highz_fwhm_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = disk_highz_fwhm_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_highz_fwhm_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = disk_highz_fwhm_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_h3(disk_highz_fwhm_BC03hr$res_h3,
        zlim = disk_highz_fwhm_BC03hr$h3_lims,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",#units = expression("residual, h"[3]), 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_BC03hr$res_h3,  yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.65,0.65), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -22, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

#h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_highz_fwhm_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = disk_highz_fwhm_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_highz_fwhm_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = disk_highz_fwhm_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_h4(disk_highz_fwhm_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",#units = expression("residual, h"[4]),
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.65,0.65), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -19, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()



# fwhm lowz blur disk EMILES ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_4/cs4_disk_velocities_lowz_fwhm_blur_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))
# velocity
plot_velocity(disk_lowz_fwhm_blur_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = disk_lowz_fwhm_blur_EMILES$vel_lims, labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_EMILES))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm_blur_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = disk_lowz_fwhm_blur_EMILES$vel_lims, labN = 2, titleshift = -5, units="",#
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm_blur_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm_blur_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-10,10),
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_EMILES$res_velocity/EMILES_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -65, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 360, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm_blur_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_blur_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm_blur_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_blur_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_velocity(disk_lowz_fwhm_blur_EMILES$res_dispersion, zlim = c(-20,20),
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_EMILES$res_dispersion/EMILES_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -90, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm_blur_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_blur_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm_blur_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#zlim = disk_lowz_fwhm_blur_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_h3(disk_lowz_fwhm_blur_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, units="",#
        #units = expression("residual, h"[3]),
        zlim = c(-0.1,0.1),#
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_EMILES$res_h3, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.2,0.2), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -75, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm_blur_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_blur_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm_blur_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_lowz_fwhm_blur_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

plot_h4(disk_lowz_fwhm_blur_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        zlim = c(-0.1,0.1),
        #units = expression("residual, h"[4]),
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_EMILES$res_h4, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.2,0.2), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -65, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()


# fwhm lowz blur disk BC03hr ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots,  "case_study_4/cs4_disk_velocities_lowz_fwhm_blur_BC03.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

# Velocity
plot_velocity(disk_lowz_fwhm_blur_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
lines(x = c(s, s+(1/lowz_scale_BC03)), y = c(s,s), lwd=2, col="black")
text((s+(1/lowz_scale_BC03))/2 + 1, 5, labels="1 kpc",adj=c(0.5,0.5), cex=1, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_lowz_fwhm_blur_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units = "",
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(disk_lowz_fwhm_blur_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_lowz_fwhm_blur_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, #units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, units = "",
              zlim = c(-50,50),
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_BC03hr$res_velocity/BC03_obs_lowz$vbin_size, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 123, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# Dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_lowz_fwhm_blur_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_blur_BC03hr$disp_lims, labN = 1, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_lowz_fwhm_blur_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_lowz_fwhm_blur_BC03hr$disp_lims, labN = 1, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_velocity(disk_lowz_fwhm_blur_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, 
              units="",#units = expression("residual, km s"^{-1}),
              zlim = c(-100, 100),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_BC03hr$res_dispersion/BC03_obs_lowz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -18, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_lowz_fwhm_blur_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = c(-0.06,0.06), labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_lowz_fwhm_blur_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = disk_lowz_fwhm_blur_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_h3(disk_lowz_fwhm_blur_BC03hr$res_h3,
        zlim = disk_lowz_fwhm_blur_BC03hr$h3_lims,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",#units = expression("residual, h"[3]), 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_BC03hr$res_h3,  yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -22, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

#h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_lowz_fwhm_blur_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.15,0.15), labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_lowz_fwhm_blur_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = disk_lowz_fwhm_blur_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

plot_h4(disk_lowz_fwhm_blur_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",#units = expression("residual, h"[4]),
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_lowz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_lowz_fwhm_blur_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.5,0.5), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -19, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()



# fwhm highz blur disk EMILES ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_4/cs4_disk_velocities_highz_fwhm_blur_EMILES.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))
# velocity
plot_velocity(disk_highz_fwhm_blur_EMILES$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
lines(x = c(s+2, s+2+(10/highz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+4+(10/highz_scale_EMILES))/2, 5, labels="10 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(15, 33, labels="Velocity Cubes",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-5, 15, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_highz_fwhm_blur_EMILES$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",#
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(15, 33, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(19, 27, labels=bquote(chi^2/DOF == .(round(median(disk_highz_fwhm_blur_EMILES$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_highz_fwhm_blur_EMILES$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, zlim = c(-10,10),
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(15, 33, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_EMILES$res_velocity/EMILES_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -65, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)
text(0, 465, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_highz_fwhm_blur_EMILES$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_blur_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_highz_fwhm_blur_EMILES$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_blur_EMILES$disp_lims, labN = 2, titleshift = -5, units="",#
                radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_velocity(disk_highz_fwhm_blur_EMILES$res_dispersion, zlim = c(-20,20),
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, units="",#units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_EMILES$res_dispersion/EMILES_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -90, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_highz_fwhm_blur_EMILES$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#disk_highz_fwhm_blur_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_highz_fwhm_blur_EMILES$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, units="",#
        zlim = c(-0.1,0.1),#zlim = disk_highz_fwhm_blur_EMILES$h3_lims, 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_h3(disk_highz_fwhm_blur_EMILES$res_h3,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, units="",#
        #units = expression("residual, h"[3]),
        zlim = c(-0.1,0.1),#
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_EMILES$res_h3, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -75, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_highz_fwhm_blur_EMILES$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_highz_fwhm_blur_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)
text(-5, 15, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_highz_fwhm_blur_EMILES$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = c(-0.1,0.1),#disk_highz_fwhm_blur_EMILES$h4_lims, 
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

plot_h4(disk_highz_fwhm_blur_EMILES$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        zlim = c(-0.1,0.1),
        #units = expression("residual, h"[4]),
        labN = 2, titleshift = -5, units="",#
        radii = disk_ellipse_data_EMILES_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_EMILES$res_h4, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.4,0.4), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -65, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()




# fwhm highz blur disk BC03hr ------------------------------------------------------------------
jpeg(filename = paste0(dir_plots, "case_study_4/cs4_disk_velocities_highz_fwhm_blur_BC03.jpeg"), 
     width = 880, height = 880, res = 100)
par(mar = c(3,2.5,2.5,0.5))

# Velocity
plot_velocity(disk_highz_fwhm_blur_BC03hr$obs_velocity, fig = c(w_array[1,], w_array[4,]), 
              zlim = c(-100,100), labN = 2, titleshift = -5, units="",
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
lines(x = c(s, s+(10/highz_scale_EMILES)), y = c(s,s), lwd=2, col="black")
text((s+(10/highz_scale_EMILES))/2 + 1, 5, labels="10 kpc",adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(12, 26, labels="Velocity Cubes", adj=c(0.5,0.5), cex=1.2, xpd=NA)
text(-4, 12, labels=expression("v"[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_velocity(disk_highz_fwhm_blur_BC03hr$ppxf_velocity, fig = c(w_array[2,], w_array[4,]), new = T, 
              zlim = c(-100,100), labN = 2, titleshift = -5, units = "",
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(12, 26, labels="Spectral Cubes",adj=c(0.5,0.5),cex=1.2, xpd=NA)
#text(15, 22, labels=bquote(chi^2/DOF == .(round(median(disk_highz_fwhm_blur_BC03hr$ppxf_chi2, na.rm=T), digits = 2))), adj=c(0.5,0.5),cex=1, xpd=NA)

plot_velocity(disk_highz_fwhm_blur_BC03hr$res_velocity,
              fig = c(w_array[3,], w_array[4,]), new = T, legend = T, #units = expression("residual, km s"^{-1}),
              labN = 2, titleshift = -5, units = "",
              zlim = c(-50,50),
              radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(12, 26, labels="Kinematic - Spectral",adj=c(0.5,0.5),cex=1.2, xpd=NA)

par(fig = c(w_array[4,], w_array[4,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_BC03hr$res_velocity/BC03_obs_highz$vbin_size, yaxt="n", 
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -20, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
abline(v=0)
text(0, 182, labels="Residual",adj=c(0.5,0.5), xpd=NA)

# Dispersion
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_dispersion(disk_highz_fwhm_blur_BC03hr$obs_dispersion, fig = c(w_array[1,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_blur_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression(paste(sigma)[LOS]*", km s"^{-1}), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_dispersion(disk_highz_fwhm_blur_BC03hr$ppxf_dispersion, fig = c(w_array[2,], w_array[3,]), new = T, 
                zlim = disk_highz_fwhm_blur_BC03hr$disp_lims, labN = 2, titleshift = -5, units="",
                radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_velocity(disk_highz_fwhm_blur_BC03hr$res_dispersion,
              fig = c(w_array[3,], w_array[3,]), new = T, legend = T, 
              units="",#units = expression("residual, km s"^{-1}),
              zlim = c(-100, 100),
              labN = 2, titleshift = -5,
              radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[3,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_BC03hr$res_dispersion/BC03_obs_highz$vbin_size, yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-1, 1), mtline = 0.5, mgp=c(1,0.25,0))
legend("topright", legend=c(expression("/"*paste(Delta)*paste(nu))), bty="n", cex=1.2)
#text(0, -18, labels=expression("residual, km s"^{-1}),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

# h3
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h3(disk_highz_fwhm_blur_BC03hr$obs_h3, fig = c(w_array[1,], w_array[2,]), new = T, 
        zlim = disk_highz_fwhm_blur_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression("h"[3]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h3(disk_highz_fwhm_blur_BC03hr$ppxf_h3, fig = c(w_array[2,], w_array[2,]), new = T, 
        zlim = disk_highz_fwhm_blur_BC03hr$h3_lims, labN = 2, titleshift = -5, units = "",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_h3(disk_highz_fwhm_blur_BC03hr$res_h3,
        zlim = disk_highz_fwhm_blur_BC03hr$h3_lims,
        fig = c(w_array[3,], w_array[2,]), new = T, legend = T, 
        units = "",#units = expression("residual, h"[3]), 
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[2,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_BC03hr$res_h3,  yaxt = "n",
        col = "#785EF0", border = NA, xlim = c(-0.65,0.65), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -22, labels=expression("residual, h"[3]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

#h4
par(mar = c(3,2.5,2.5,0.5), cex=1)
plot_h4(disk_highz_fwhm_blur_BC03hr$obs_h4, fig = c(w_array[1,], w_array[1,]), new = T, 
        zlim = disk_highz_fwhm_blur_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)
text(-4, 12, labels=expression("h"[4]), adj=c(0.5,0.5), cex=1.4, xpd=NA, srt=90)

plot_h4(disk_highz_fwhm_blur_BC03hr$ppxf_h4, fig = c(w_array[2,], w_array[1,]), new = T, 
        zlim = disk_highz_fwhm_blur_BC03hr$h4_lims, labN = 2, titleshift = -5,units="",
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

plot_h4(disk_highz_fwhm_blur_BC03hr$res_h4,
        fig = c(w_array[3,], w_array[1,]), new = T, legend = T, 
        units="",#units = expression("residual, h"[4]),
        labN = 2, titleshift = -5,
        radii = disk_ellipse_data_BC03_highz, cex=1.2)

par(fig = c(w_array[4,],w_array[1,]), new=T, cex=1.2)
maghist(disk_highz_fwhm_blur_BC03hr$res_h4, yaxt="n",
        col = "#785EF0", border = NA, xlim = c(-0.65,0.65), mtline = 0.5, mgp=c(1,0.25,0))
#text(0, -19, labels=expression("residual, h"[4]),adj=c(0.5,0.5),cex=1, xpd=NA)
abline(v=0)

dev.off()



# histogram summary cs3 ------------------------------------------------------------
plot_dist = function(data_for_hist, lty, col){
  h = maghist(data_for_hist, breaks = seq(-1, 1, length.out=50), plot = F)
  lines(h$mids, h$counts, lwd=2, lty=lty, col=col)
}

w_array = plot_array(plot_size = 0.55, n_plots = 2)
h_array = plot_array(plot_size = 0.54, n_plots = 2)

jpeg(filename = paste0(dir_plots, "case_study_3/cs3_histograms.jpeg"), 
     width = 610, height = 650, res = 100)
par(mar = c(4.5,3.5,1.5,0.5))

par(fig = c(w_array[1,], h_array[2,]))
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,200),
        xlab=expression("residual v"["LOS"]*"/"*paste(Delta)*paste(nu)),
        ylab="Pixels")
plot_dist(disk_lowz_fwhm_EMILES$res_velocity/EMILES_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_BC03hr$res_velocity/BC03_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_EMILES$res_velocity/EMILES_obs_highz$vbin_size, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_BC03hr$res_velocity/BC03_obs_highz$vbin_size, 
          col = "#40B0A6", lty=3)
abline(v=0)
legend("topleft", legend = c("z = 0.0144", "z = 0.3"),
       col = c("#E1BE6A", "#40B0A6"), pch=15, bty = "n")


par(fig = c(w_array[2,], h_array[2,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,200),
        xlab=expression("residual "*paste(sigma)["LOS"]*"/"*paste(Delta)*paste(nu)),
        yaxt="n")
plot_dist(disk_lowz_fwhm_EMILES$res_dispersion/EMILES_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_BC03hr$res_dispersion/BC03_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_EMILES$res_dispersion/EMILES_obs_highz$vbin_size, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_BC03hr$res_dispersion/BC03_obs_highz$vbin_size, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

legend("topright", legend = c("BC03hr", "EMILES"),
       col = c("black", "black"), lty=c(3,1), lwd=3, bty = "n")

par(fig = c(w_array[1,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,400),
        xlab=expression("residual h"["3"]),
        ylab="Pixels")
plot_dist(disk_lowz_fwhm_EMILES$res_h3, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_BC03hr$res_h3, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_EMILES$res_h3, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_BC03hr$res_h3, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

par(fig = c(w_array[2,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,400),
        xlab=expression("residual h"["4"]),
        yaxt="n")
plot_dist(disk_lowz_fwhm_EMILES$res_h4, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_BC03hr$res_h4, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_EMILES$res_h4, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_BC03hr$res_h4, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

dev.off()

# histogram summary cs4 ------------------------------------------------------------
plot_dist = function(data_for_hist, lty, col){
  h = maghist(data_for_hist, breaks = seq(-1, 1, length.out=50), plot = F)
  lines(h$mids, h$counts, lwd=2, lty=lty, col=col)
}

w_array = plot_array(plot_size = 0.55, n_plots = 2)
h_array = plot_array(plot_size = 0.54, n_plots = 2)

jpeg(filename = paste0(dir_plots, "case_study_4/cs4_histograms.jpeg"), 
     width = 610, height = 650, res = 100)
par(mar = c(4.5,3.5,1.5,0.5))

par(fig = c(w_array[1,], h_array[2,]))
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,230),
        xlab=expression("residual v"["LOS"]*"/"*paste(Delta)*paste(nu)),
        ylab="Pixels")
plot_dist(disk_lowz_fwhm_blur_EMILES$res_velocity/EMILES_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_blur_BC03hr$res_velocity/BC03_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_blur_EMILES$res_velocity/EMILES_obs_highz$vbin_size, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_blur_BC03hr$res_velocity/BC03_obs_highz$vbin_size, 
          col = "#40B0A6", lty=3)
abline(v=0)
legend("topleft", legend = c("z = 0.0144", "z = 0.3"),
       col = c("#E1BE6A", "#40B0A6"), pch=15, bty = "n")

par(fig = c(w_array[2,], h_array[2,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,230),
        xlab=expression("residual "*paste(sigma)["LOS"]*"/"*paste(Delta)*paste(nu)),
        yaxt="n")
plot_dist(disk_lowz_fwhm_blur_EMILES$res_dispersion/EMILES_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_blur_BC03hr$res_dispersion/BC03_obs_lowz$vbin_size, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_blur_EMILES$res_dispersion/EMILES_obs_highz$vbin_size, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_blur_BC03hr$res_dispersion/BC03_obs_highz$vbin_size, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

legend("topright", legend = c("BC03hr", "EMILES"),
       col = c("black", "black"), lty=c(3,1), lwd=3, bty = "n")

par(fig = c(w_array[1,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,400),
        xlab=expression("residual h"["3"]),
        ylab="Pixels")
plot_dist(disk_lowz_fwhm_blur_EMILES$res_h3, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_blur_BC03hr$res_h3, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_blur_EMILES$res_h3, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_blur_BC03hr$res_h3, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

par(fig = c(w_array[2,], h_array[1,]), new=T)
magplot(0, type="n", xlim = c(-1,1), ylim = c(0,400),
        xlab=expression("residual h"["4"]),
        yaxt="n")
plot_dist(disk_lowz_fwhm_blur_EMILES$res_h4, 
          col = "#E1BE6A", lty=1)
plot_dist(disk_lowz_fwhm_blur_BC03hr$res_h4, 
          col = "#E1BE6A", lty=3)
plot_dist(disk_highz_fwhm_blur_EMILES$res_h4, 
          col = "#40B0A6", lty=1)
plot_dist(disk_highz_fwhm_blur_BC03hr$res_h4, 
          col = "#40B0A6", lty=3)
abline(v=0, lwd=2)

dev.off()
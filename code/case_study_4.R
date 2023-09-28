# Comparison of spectral and kinematic cubes for Case Study 4
# Kate Harborne 
#
# Running this file will generate the data for constructing the plots in Case 
# Study 4 (section 3.1 test 4)
#
################################################################################
#
# Installation of packages below may be necessary if you do not have them 
# already installed. Uncomment the lines below to install the necessary modules.

# install.packages("devtools")
# devtools::install_github("kateharborne/SimSpin@v2.6.0")
# install.packages("magicaxis")
#
################################################################################

# Loading packges -------------------------------------------------------------- 

library(SimSpin)
packageVersion("SimSpin")
library(magicaxis)

# Reference directories --------------------------------------------------------

ss_file_dir     = "data/simspin_files/"
ss_datacube_dir_1 = "data/cubes/fwhm_lowz_blur/"
ss_datacube_dir_2 = "data/cubes/fwhm_highz_blur/"

# Initialising telescope and observing strategy --------------------------------

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
                                   blur=T, fwhm = 1, 
                                   psf="Gaussian")

strategy_two  = observing_strategy(dist_z = 0.3, 
                                   inc_deg = 60, 
                                   blur=T, fwhm = 2.8, 
                                   psf="Moffat")


# SPECTRAL CUBES ===============================================================
# Case Study 3: (1a) Observations at low redshift EMILES (DISK) ----

disk_age05_Z024_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                          telescope = SAMI, 
                                          observing_strategy = strategy_one,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir_1, "/disk_EMILES/disk_EMILES_spectral_fwhm_lowz_blur.FITS"),
                                          object_name = "disk_age05_Z004_EMILES", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_spectral, paste0(ss_datacube_dir_1, "/disk_EMILES/disk_EMILES_spectral_fwhm_lowz_blur.Rdata"))
remove(disk_age05_Z024_spectral)

# Case Study 3: (1b) Observations at high redshift EMILES (DISK) ----

disk_age05_Z024_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                          telescope = SAMI, 
                                          observing_strategy = strategy_two,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir_2, "/disk_EMILES/disk_EMILES_spectral_fwhm_highz_blur.FITS"),
                                          object_name = "disk_age05_Z004_EMILES", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_spectral, paste0(ss_datacube_dir_2, "/disk_EMILES/disk_EMILES_spectral_fwhm_highz_blur.Rdata"))
remove(disk_age05_Z024_spectral)

# Case Study 3: (2a) Observations at low redshift (disk) with BC03 ----

disk_age05_Z024_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                          telescope = SAMI_bc, 
                                          observing_strategy = strategy_one,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir_1, "/disk_BC03hr/disk_BC03hr_spectral_fwhm_lowz_blur.FITS"),
                                          object_name = "disk_age05_Z004_BC03hr", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_spectral, paste0(ss_datacube_dir_1, "/disk_BC03hr/disk_BC03hr_spectral_fwhm_lowz_blur.Rdata"))
remove(disk_age05_Z024_spectral)

# Case Study 3: (2b) Observations at high redshift with BC03 (DISK) ----

disk_age05_Z024_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                          telescope = SAMI_bc, 
                                          observing_strategy = strategy_two,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir_2, "/disk_BC03hr/disk_BC03hr_spectral_fwhm_highz_blur.FITS"),
                                          object_name = "disk_age05_Z004_BC03hr", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_spectral, paste0(ss_datacube_dir_2, "/disk_BC03hr/disk_BC03hr_spectral_fwhm_highz_blur.Rdata"))
remove(disk_age05_Z024_spectral)

# KINEMATIC CUBES ===============================================================
# Case Study 3: (1a) Observations at low redshift EMILES (DISK) ----

disk_age05_Z024_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                           telescope = SAMI, 
                                           observing_strategy = strategy_one,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir_1, "/disk_EMILES/disk_EMILES_kinematic_fwhm_lowz_blur.FITS"),
                                           object_name = "disk_age05_Z004_EMILES", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_kinematic, paste0(ss_datacube_dir_1, "/disk_EMILES/disk_EMILES_kinematic_fwhm_lowz_blur.Rdata"))
remove(disk_age05_Z024_kinematic)

# Case Study 3: (1b) Observations at high redshift EMILES (DISK) ----

disk_age05_Z024_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                           telescope = SAMI, 
                                           observing_strategy = strategy_two,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir_2, "/disk_EMILES/disk_EMILES_kinematic_fwhm_highz_blur.FITS"),
                                           object_name = "disk_age05_Z004_EMILES", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_kinematic, paste0(ss_datacube_dir_2, "/disk_EMILES/disk_EMILES_kinematic_fwhm_highz_blur.Rdata"))
remove(disk_age05_Z024_kinematic)

# Case Study 3: (2a) Observations at low redshift (disk) with BC03 ----

disk_age05_Z024_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                           telescope = SAMI_bc, 
                                           observing_strategy = strategy_one,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir_1, "/disk_BC03hr/disk_BC03hr_kinematic_fwhm_lowz_blur.FITS"),
                                           object_name = "disk_age05_Z004_BC03hr", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_kinematic, paste0(ss_datacube_dir_1, "/disk_BC03hr/disk_BC03hr_kinematic_fwhm_lowz_blur.Rdata"))
remove(disk_age05_Z024_kinematic)

# Case Study 3: (2b) Observations at high redshift with BC03 (DISK) ----

disk_age05_Z024_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                           telescope = SAMI_bc, 
                                           observing_strategy = strategy_two,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir_2, "/disk_BC03hr/disk_BC03hr_kinematic_fwhm_highz_blur.FITS"),
                                           object_name = "disk_age05_Z004_BC03hr", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z024_kinematic, paste0(ss_datacube_dir_2, "/disk_BC03hr/disk_BC03hr_kinematic_fwhm_highz_blur.Rdata"))
remove(disk_age05_Z024_kinematic)

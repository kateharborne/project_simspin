# Comparison of spectral and kinematic cubes for Case Study 1
# Kate Harborne 
#
# Running this file will generate the data for constructing the plots in Case 
# Study 1 (section 3.1 test 1)
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
ss_datacube_dir = "data/cubes/fwhm0_lowz/"

# Initialising telescope and observing strategy --------------------------------

SAMI = telescope(type = "IFU", 
                 signal_to_noise = 30, 
                 lsf_fwhm = 0, 
                 fov = 15, 
                 wave_range = c(3750,5750))

SAMI_bc = telescope(type = "IFU", 
                    signal_to_noise = 30, 
                    lsf_fwhm = 0, 
                    spatial_res = 0.7, 
                    wave_res = 3,
                    fov = 17, 
                    wave_range = c(3750,5750),
                    aperture_shape = "hexagonal")

strategy_one  = observing_strategy(dist_kpc_per_arcsec = 0.3, 
                                   inc_deg = 60, 
                                   blur=F)

# Generating the SimSpin files necessary for these tests -----------------------

make_simspin_file(filename = "data/galaxy_cutouts/disk.hdf5",
                  disk_age = 5, disk_Z = 0.004,   
                  output = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"), template = "EMILES",
                  overwrite = T)

make_simspin_file(filename = "data/galaxy_cutouts/disk.hdf5",
                  disk_age = 5, disk_Z = 0.004,   
                  output = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"), template = "BC03hr",
                  overwrite = T)

make_simspin_file(filename = "data/galaxy_cutouts/bulge.hdf5",
                  bulge_age = 5, bulge_Z = 0.004, 
                  output = paste0(ss_file_dir, "/bulge_age05_Z004_EMILES.Rdata"), template = "EMILES",
                  overwrite = T)

make_simspin_file(filename = "data/galaxy_cutouts/bulge.hdf5",
                  bulge_age = 5, bulge_Z = 0.004, 
                  output = paste0(ss_file_dir, "/bulge_age05_Z004_BC03hr.Rdata"), template = "BC03hr",
                  overwrite = T)

make_simspin_file(filename = "data/galaxy_cutouts/bulge.hdf5",
                  bulge_age = 10, bulge_Z = 0.001, 
                  output = paste0(ss_file_dir, "/bulge_age10_Z001_EMILES.Rdata"), template = "EMILES",
                  overwrite = T)

make_simspin_file(filename = "data/galaxy_cutouts/bulge.hdf5",
                  bulge_age = 10, bulge_Z = 0.001, 
                  output = paste0(ss_file_dir, "/bulge_age10_Z001_BC03hr.Rdata"), template = "BC03hr",
                  overwrite = T)

# SPECTRAL CUBES ===============================================================
# Test 1: (1a) Low redshift DISK built with EMILES at fwhm 0 -------------------

disk_age05_Z004_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                          telescope = SAMI, 
                                          observing_strategy = strategy_one,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir, "/disk_EMILES/disk_EMILES_spectral_fwhm0_lowz.FITS"),
                                          object_name = "disk_age05_Z004", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z004_spectral, paste0(ss_datacube_dir, "/disk_EMILES/disk_EMILES_spectral_fwhm0_lowz.Rdata"))
remove(disk_age05_Z004_spectral)

# Test 1: (2a) Low redshift BULGE built with EMILES at fwhm 0 ------------------

bulge_age05_Z004_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age05_Z004_EMILES.Rdata"),
                                           telescope = SAMI, 
                                           observing_strategy = strategy_one,
                                           method = "spectral",
                                           write_fits = T, cores = 2, split_save = F,
                                           output_location = paste0(ss_datacube_dir, "/bulge_EMILES/bulge_EMILES_spectral_fwhm0_lowz.FITS"),
                                           object_name = "bulge_age05_Z004", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age05_Z004_spectral, paste0(ss_datacube_dir, "/bulge_EMILES/bulge_EMILES_spectral_fwhm0_lowz.Rdata"))
remove(bulge_age05_Z004_spectral)

# Test 1: (3a) Low redshift OLD BULGE with EMILES at z = 0 ----

bulge_age10_Z001_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age10_Z001_EMILES.Rdata"),
                                           telescope = SAMI, 
                                           observing_strategy = strategy_one,
                                           method = "spectral",
                                           write_fits = T, cores = 2, split_save = F,
                                           output_location = paste0(ss_datacube_dir, "/old_bulge_EMILES/old_bulge_EMILES_spectral_fwhm0_lowz.FITS"),
                                           object_name = "bulge_age10_Z001", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age10_Z001_spectral, paste0(ss_datacube_dir, "/old_bulge_EMILES/old_bulge_EMILES_spectral_fwhm0_lowz.Rdata"))
remove(bulge_age10_Z001_spectral)

# Test 1: (1b) Low redshift DISK built with BC03hr at fwhm 0 -------------------

disk_age05_Z004_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                          telescope = SAMI_bc, 
                                          observing_strategy = strategy_one,
                                          method = "spectral",
                                          write_fits = T, cores = 2, split_save = F,
                                          output_location = paste0(ss_datacube_dir, "/disk_BC03hr/disk_BC03hr_spectral_fwhm0_lowz.FITS"),
                                          object_name = "disk_age05_Z004", 
                                          telescope_name = "SAMI", 
                                          observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z004_spectral, paste0(ss_datacube_dir, "/disk_BC03hr/disk_BC03hr_spectral_fwhm0_lowz.Rdata"))
remove(disk_age05_Z004_spectral)

# Test 1: (2b) Low redshift BULGE built with BC03hr at fwhm 0 ------------------

bulge_age05_Z004_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age05_Z004_BC03hr.Rdata"),
                                           telescope = SAMI_bc, 
                                           observing_strategy = strategy_one,
                                           method = "spectral",
                                           write_fits = T, cores = 2, split_save = F,
                                           output_location = paste0(ss_datacube_dir, "/bulge_BC03hr/bulge_BC03hr_spectral_fwhm0_lowz.FITS"),
                                           object_name = "bulge_age05_Z004", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age05_Z004_spectral, paste0(ss_datacube_dir, "/bulge_BC03hr/bulge_BC03hr_spectral_fwhm0_lowz.Rdata"))
remove(bulge_age05_Z004_spectral)

# Test 1: (3b) Low redshift OLD BULGE with BC03hr at z = 0 ---------------------

bulge_age10_Z001_spectral = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age10_Z001_BC03hr.Rdata"),
                                           telescope = SAMI_bc, 
                                           observing_strategy = strategy_one,
                                           method = "spectral",
                                           write_fits = T, cores = 2, split_save = F,
                                           output_location = paste0(ss_datacube_dir, "/old_bulge_BC03hr/old_bulge_BC03hr_spectral_fwhm0_lowz.FITS"),
                                           object_name = "bulge_age10_Z001", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age10_Z001_spectral, paste0(ss_datacube_dir, "/old_bulge_BC03hr/old_bulge_BC03hr_spectral_fwhm0_lowz.Rdata"))
remove(bulge_age10_Z001_spectral)

# KINEMATIC CUBES ==============================================================

# Test 1: (1a) Low redshift DISK built with EMILES at fwhm 0 -------------------

disk_age05_Z004_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_EMILES.Rdata"),
                                           telescope = SAMI, 
                                           observing_strategy = strategy_one,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir, "/disk_EMILES/disk_EMILES_kinematic_fwhm0_lowz.FITS"),
                                           object_name = "disk_age05_Z004", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z004_kinematic, paste0(ss_datacube_dir,"/disk_EMILES/disk_EMILES_kinematic_fwhm0_lowz.Rdata"))
remove(disk_age05_Z004_kinematic)

# Test 1: (2a) Low redshift BULGE built with EMILES at fwhm 0 ------------------

bulge_age05_Z004_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age05_Z004_EMILES.Rdata"),
                                            telescope = SAMI, 
                                            observing_strategy = strategy_one,
                                            method = "velocity",
                                            write_fits = T, cores = 2, 
                                            output_location = paste0(ss_datacube_dir, "/bulge_EMILES/bulge_EMILES_kinematic_fwhm0_lowz.FITS"),
                                            object_name = "bulge_age05_Z004", 
                                            telescope_name = "SAMI", 
                                            observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age05_Z004_kinematic, paste0(ss_datacube_dir,"/bulge_EMILES/bulge_EMILES_kinematic_fwhm0_lowz.Rdata"))
remove(bulge_age05_Z004_kinematic)

# Test 1: (3a) Low redshift OLD BULGE with EMILES at z = 0 ---------------------

bulge_age10_Z001_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age10_Z001_EMILES.Rdata"),
                                            telescope = SAMI, 
                                            observing_strategy = strategy_one,
                                            method = "velocity",
                                            write_fits = T, cores = 2, 
                                            output_location = paste0(ss_datacube_dir, "/old_bulge_EMILES/old_bulge_EMILES_kinematic_fwhm0_lowz.FITS"),
                                            object_name = "bulge_age10_Z001", 
                                            telescope_name = "SAMI", 
                                            observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age10_Z001_kinematic, paste0(ss_datacube_dir,"/old_bulge_EMILES/old_bulge_EMILES_kinematic_fwhm0_lowz.Rdata"))
remove(bulge_age10_Z001_kinematic)

# Test 1: (1b) Low redshift DISK built with BC03hr at fwhm 0 -------------------

disk_age05_Z004_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/disk_age05_Z004_BC03hr.Rdata"),
                                           telescope = SAMI_bc, 
                                           observing_strategy = strategy_one,
                                           method = "velocity",
                                           write_fits = T, cores = 2, 
                                           output_location = paste0(ss_datacube_dir, "/disk_BC03hr/disk_BC03hr_kinematic_fwhm0_lowz.FITS"),
                                           object_name = "disk_age05_Z004", 
                                           telescope_name = "SAMI", 
                                           observer_name = "K E Harborne", verbose = F)

saveRDS(disk_age05_Z004_kinematic, paste0(ss_datacube_dir,"/disk_BC03hr/disk_BC03hr_kinematic_fwhm0_lowz.Rdata"))
remove(disk_age05_Z004_kinematic)


# Test 1: (2b) Low redshift BULGE built with BC03hr at fwhm 0 ------------------

bulge_age05_Z004_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age05_Z004_BC03hr.Rdata"),
                                            telescope = SAMI_bc, 
                                            observing_strategy = strategy_one,
                                            method = "velocity",
                                            write_fits = T, cores = 2, 
                                            output_location = paste0(ss_datacube_dir, "/bulge_BC03hr/bulge_BC03hr_kinematic_fwhm0_lowz.FITS"),
                                            object_name = "bulge_age05_Z004", 
                                            telescope_name = "SAMI", 
                                            observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age05_Z004_kinematic, paste0(ss_datacube_dir,"/bulge_BC03hr/bulge_BC03hr_kinematic_fwhm0_lowz.Rdata"))
remove(bulge_age05_Z004_kinematic)


# Test 1: (3b) Low redshift OLD BULGE with BC03hr at z = 0 ---------------------

bulge_age10_Z001_kinematic = build_datacube(simspin_file = paste0(ss_file_dir, "/bulge_age10_Z001_BC03hr.Rdata"),
                                            telescope = SAMI_bc, 
                                            observing_strategy = strategy_one,
                                            method = "velocity",
                                            write_fits = T, cores = 2, 
                                            output_location = paste0(ss_datacube_dir, "/old_bulge_BC03hr/old_bulge_BC03hr_kinematic_fwhm0_lowz.FITS"),
                                            object_name = "bulge_age10_Z001", 
                                            telescope_name = "SAMI", 
                                            observer_name = "K E Harborne", verbose = F)

saveRDS(bulge_age10_Z001_kinematic, paste0(ss_datacube_dir,"/old_bulge_BC03hr/old_bulge_BC03hr_kinematic_fwhm0_lowz.Rdata"))
remove(bulge_age10_Z001_kinematic)


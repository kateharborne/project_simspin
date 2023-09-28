library(reticulate)
np = import("numpy", convert=T)  

summarise_results = function(fname, folder, plot = T){
  # Reading in the kinematic observation
  simspin_kincube = readRDS(paste0(simspin_out, folder, stringr::str_replace(fname, "spectral", "kinematic"),".Rdata"))
  
  # Determine the aperture
  aperture = matrix(data = simspin_kincube$observation$aperture_region, 
                    nrow = dim(simspin_kincube$raw_images$particle_image)[1], 
                    ncol = dim(simspin_kincube$raw_images$particle_image)[1])
  aperture[aperture == 0] = NA
  
  if (plot){
    par(fig = c(0,0.5,0,1))
    maghist(simspin_kincube$observed_images$residuals)
    plot_age(simspin_kincube$observed_images$residuals, fig = c(0.5,1,0,1),
             new=T, labN = 2, titleshift = -5, units = "LOSVD Residuals")
  }
  
  #masking pixels with >5% error on fitted losvd
  aperture[simspin_kincube$observed_images$residuals > 0.05] = NA
  
  summary_kin = 
    list(ppxf_velocity   = t(np$load(paste0(ppxf_out, folder, fname, "_ppxf_velocity.npy"))), 
         ppxf_dispersion = t(np$load(paste0(ppxf_out, folder, fname, "_ppxf_dispersion.npy"))), 
         ppxf_h3  = t(np$load(paste0(ppxf_out, folder, fname, "_ppxf_h3.npy"))), 
         ppxf_h4 = t(np$load(paste0(ppxf_out, folder, fname,"_ppxf_h4.npy"))), 
         ppxf_chi2 = t(np$load(paste0(ppxf_out, folder, fname, "_ppxf_chi2.npy"))))
         
  summary_kin$obs_velocity = (simspin_kincube$observed_images$velocity_image - median(simspin_kincube$observed_images$velocity_image, na.rm = T)) * aperture 
  summary_kin$obs_dispersion = simspin_kincube$observed_images$dispersion_image * aperture
  summary_kin$obs_h3 = (simspin_kincube$observed_images$h3_image - median(simspin_kincube$observed_images$h3_image, na.rm = T)) * aperture
  summary_kin$obs_h4 = (simspin_kincube$observed_images$h4_image - median(simspin_kincube$observed_images$h4_image, na.rm = T)) * aperture
  summary_kin$ppxf_velocity = (summary_kin$ppxf_velocity - median(summary_kin$ppxf_velocity, na.rm = T)) * aperture
  summary_kin$ppxf_dispersion = summary_kin$ppxf_dispersion * aperture
  summary_kin$ppxf_h3 = (summary_kin$ppxf_h3 - median(summary_kin$ppxf_h3, na.rm = T)) * aperture
  summary_kin$ppxf_h4 = (summary_kin$ppxf_h4 - median(summary_kin$ppxf_h4, na.rm = T)) * aperture
  summary_kin$num_part = simspin_kincube$raw_images$particle_image * aperture
  summary_kin$radius = get_radius(simspin_kincube, aperture)
  summary_kin$res_velocity = summary_kin$obs_velocity - summary_kin$ppxf_velocity
  summary_kin$res_dispersion = summary_kin$obs_dispersion - summary_kin$ppxf_dispersion
  summary_kin$res_h3 = summary_kin$obs_h3 - summary_kin$ppxf_h3
  summary_kin$res_h4 = summary_kin$obs_h4 - summary_kin$ppxf_h4
  
  # Computing summaries of distribution as function of RADIUS
  summary_kin$r_bins = cut(summary_kin$radius, 
                           seq(0, max(summary_kin$radius, na.rm=T)+0.5, by = 0.5),
                           labels = F)
  summary_kin$r_bin_labels = seq(0, max(summary_kin$radius, na.rm=T)+0.5, by = 0.5)
  
  summary_kin$res_r_velocity_med = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_dispersion_med = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h3_med = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h4_med = numeric(length(summary_kin$r_bin_labels))
  
  summary_kin$res_r_velocity_16 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_dispersion_16 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h3_16 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h4_16 = numeric(length(summary_kin$r_bin_labels))
  
  summary_kin$res_r_velocity_84 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_dispersion_84 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h3_84 = numeric(length(summary_kin$r_bin_labels))
  summary_kin$res_r_h4_84 = numeric(length(summary_kin$r_bin_labels))
  
  for (each in 1:length(summary_kin$r_bin_labels)){
    summary_kin$res_r_velocity_med[each] = median(summary_kin$res_velocity[summary_kin$r_bins == each], na.rm = T)
    summary_kin$res_r_dispersion_med[each] = median(summary_kin$res_dispersion[summary_kin$r_bins == each], na.rm = T)
    summary_kin$res_r_h3_med[each] = median(summary_kin$res_h3[summary_kin$r_bins == each], na.rm = T)
    summary_kin$res_r_h4_med[each] = median(summary_kin$res_h4[summary_kin$r_bins == each], na.rm = T)
    
    summary_kin$res_r_velocity_16[each] = quantile(summary_kin$res_velocity[summary_kin$r_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_r_dispersion_16[each] = quantile(summary_kin$res_dispersion[summary_kin$r_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_r_h3_16[each] = quantile(summary_kin$res_h3[summary_kin$r_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_r_h4_16[each] = quantile(summary_kin$res_h4[summary_kin$r_bins == each], probs = c(0.16), na.rm = T)
    
    summary_kin$res_r_velocity_84[each] = quantile(summary_kin$res_velocity[summary_kin$r_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_r_dispersion_84[each] = quantile(summary_kin$res_dispersion[summary_kin$r_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_r_h3_84[each] = quantile(summary_kin$res_h3[summary_kin$r_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_r_h4_84[each] = quantile(summary_kin$res_h4[summary_kin$r_bins == each], probs = c(0.84), na.rm = T)
  }
  
  # Computing summaries of distributions as a function of PARTICLE/PIXEL
  summary_kin$n_bins = cut(summary_kin$num_part, 
                           seq(min(summary_kin$num_part, na.rm=T)-50, max(summary_kin$num_part, na.rm=T)+50, by = 50),
                           labels = F)
  summary_kin$n_bin_labels = seq(min(summary_kin$num_part, na.rm=T)-50, max(summary_kin$num_part, na.rm=T)+50, by = 50)
  
  summary_kin$res_n_velocity_med = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_dispersion_med = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h3_med = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h4_med = numeric(length(summary_kin$n_bin_labels))
  
  summary_kin$res_n_velocity_16 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_dispersion_16 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h3_16 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h4_16 = numeric(length(summary_kin$n_bin_labels))
  
  summary_kin$res_n_velocity_84 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_dispersion_84 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h3_84 = numeric(length(summary_kin$n_bin_labels))
  summary_kin$res_n_h4_84 = numeric(length(summary_kin$n_bin_labels))
  
  for (each in 1:length(summary_kin$n_bin_labels)){
    summary_kin$res_n_velocity_med[each] = median(summary_kin$res_velocity[summary_kin$n_bins == each], na.rm = T)
    summary_kin$res_n_dispersion_med[each] = median(summary_kin$res_dispersion[summary_kin$n_bins == each], na.rm = T)
    summary_kin$res_n_h3_med[each] = median(summary_kin$res_h3[summary_kin$n_bins == each], na.rm = T)
    summary_kin$res_n_h4_med[each] = median(summary_kin$res_h4[summary_kin$n_bins == each], na.rm = T)
    
    summary_kin$res_n_velocity_16[each] = quantile(summary_kin$res_velocity[summary_kin$n_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_n_dispersion_16[each] = quantile(summary_kin$res_dispersion[summary_kin$n_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_n_h3_16[each] = quantile(summary_kin$res_h3[summary_kin$n_bins == each], probs = c(0.16), na.rm = T)
    summary_kin$res_n_h4_16[each] = quantile(summary_kin$res_h4[summary_kin$n_bins == each], probs = c(0.16), na.rm = T)
    
    summary_kin$res_n_velocity_84[each] = quantile(summary_kin$res_velocity[summary_kin$n_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_n_dispersion_84[each] = quantile(summary_kin$res_dispersion[summary_kin$n_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_n_h3_84[each] = quantile(summary_kin$res_h3[summary_kin$n_bins == each], probs = c(0.84), na.rm = T)
    summary_kin$res_n_h4_84[each] = quantile(summary_kin$res_h4[summary_kin$n_bins == each], probs = c(0.84), na.rm = T)
  }
  
  # Computing the limits for each plot
  summary_kin$vel_lims = c(-max(c(abs(summary_kin$obs_velocity), abs(summary_kin$ppxf_velocity)), na.rm=T),
                           max(c(abs(summary_kin$obs_velocity), abs(summary_kin$ppxf_velocity)), na.rm=T))
  summary_kin$disp_lims = c(min(c(summary_kin$obs_dispersion, summary_kin$ppxf_dispersion), na.rm=T),
                            max(c(summary_kin$obs_dispersion, summary_kin$ppxf_dispersion), na.rm=T))
  summary_kin$h3_lims = c(-max(c(abs(summary_kin$obs_h3), abs(summary_kin$ppxf_h3)), na.rm=T),
                          max(c(abs(summary_kin$obs_h3), abs(summary_kin$ppxf_h3)), na.rm=T))
  summary_kin$h4_lims = c(-max(c(abs(summary_kin$obs_h4), abs(summary_kin$ppxf_h4)), na.rm=T),
                          max(c(abs(summary_kin$obs_h4), abs(summary_kin$ppxf_h4)), na.rm=T))
  
  return(summary_kin)
}

get_radius = function(cube, aperture){
  sbin = cube$observation$sbin
  xcoord = (cube$observation$sbin_seq + (cube$observation$sbin_size/2))[1:sbin]
  sbin_grid = list("x" = array(data = rep(xcoord, sbin), dim = c(sbin,sbin)),
                   "y" = array(data = rep(xcoord, each=sbin), dim = c(sbin,sbin)))
  radius = sqrt(sbin_grid$x^2 + sbin_grid$y^2)
  return(radius*aperture)
}
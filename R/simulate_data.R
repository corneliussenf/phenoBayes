# Simulate data

simulate_data <- function(t = 1:365, # DOY
                          a = 0.15, # Green-up
                          as = 0.005, # Pixel variation in green-up
                          b = 100, # Start of season (in DOY)
                          bs = 8, # Pixel variation in start of season (DOY)
                          c = 0.15, # Senessence
                          cs = 0.005, # Pixel variation in end of season
                          d = 270, # End of season (in DOY)
                          ds = 8, # Pixel variation in end of season (in DOY)
                          vi_0 = 0, # Minimum of VI
                          vi_0_s = 0.01, # Pixel variation in minimum of VI
                          vi_delta = 0.5, # Magnitude of VI
                          vi_delta_s = 0.01, # Pixel variation in magnitude of VI
                          l = 0.001, # Greendown parameter
                          ls = 0.0001, # Pixel variation in greendown parameter
                          error = 0.05, # Random error term
                          N = 10, # Average number of observations per year
                          S = 20, # Variation in number of observations per year
                          Y = 30, # Number of years
                          M = 25, # Number of time series
                          var_sos = 8, # Annual variation in SOS
                          var_eos = 0, # Annual variation in EOS
                          var_gup = 0, # Annual variation in GUP
                          var_gse = 0, # Annual variation in GSE
                          var_vi_0 = 0, # Annual variation in VI0
                          var_vi_delta = 0, # Annual variation in VID
                          var_gdw = 0) { # Annual variation in GDW

  # Double logistic function
  dl <- function(t, vi_0, vi_delta, a, b, c, d, l){
    l1 <- (1 / (1 + exp(- a * (t - b))))
    l2 <- (1 / (1 + exp(- c * (t - d))))
    return(vi_0 + (vi_delta - l * t) * (l1 - l2))
  }

  # Define color set
  getPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  cols <- getPalette(Y)

  # Create parameters
  a_var_y <- rnorm(Y, 0, var_gup)
  b_var_y <- rnorm(Y, 0, var_sos)
  c_var_y <- rnorm(Y, 0, var_gse)
  d_var_y <- rnorm(Y, 0, var_eos)
  vi_0_var_y <- rnorm(Y, 0, var_vi_0)
  vi_delta_var_y <- rnorm(Y, 0, var_vi_delta)
  l_var_y <- rnorm(Y, 0, var_gdw)

  dat_pred <- vector("list", M)
  dat <- vector("list", M)

  a_var <- vector("numeric", M)
  b_var <- vector("numeric", M)
  c_var <- vector("numeric", M)
  d_var <- vector("numeric", M)
  vi_0_var <- vector("numeric", M)
  vi_delta_var <- vector("numeric", M)
  l_var <- vector("numeric", M)

  n_pixels <- matrix(NA, ncol = M, nrow = Y)

  # Simulate spatial variability
  for(k in 1:M){

    n <- round(rnorm(Y, N, S), 0)
    n[n<0] <- 0
    n_pixels[, k] <- n

    a_var[k] <- rnorm(1, a, as)
    b_var[k] <- rnorm(1, b, bs)
    c_var[k] <- rnorm(1, c, cs)
    d_var[k] <- rnorm(1, d, ds)
    vi_0_var[k] <- rnorm(1, vi_0, vi_0_s)
    vi_delta_var[k] <- rnorm(1, vi_delta, vi_delta_s)
    l_var[k] <- rnorm(1, l, ls)

    vi <- vector(mode = "list", Y)

    # Simulate temporal variability
    for (i in 1:Y) {
      vi_raw <- dl(t = t,
                    vi_0 = vi_0_var[k] + vi_0_var_y[i],
                    vi_delta = vi_delta_var[k] + vi_delta_var_y[i],
                    a = a_var[k] + a_var_y[i],
                    b = b_var[k] + b_var_y[i],
                    c = c_var[k] + c_var_y[i],
                    d = d_var[k] + d_var_y[i],
                    l = l_var[k] + l_var_y[i])
      vi_error <- vi_raw + rnorm(length(vi_raw), 0, error)

      vi_error_select <- vi_error
      vi_error_select[sample(1:length(vi_error), ifelse(length(vi_error) - n[i] < 0, 0, length(vi_error) - n[i]))] <- NA
      vi_error_select <- data.frame(vi = vi_error_select, y = i, doy = t)
      vi[[i]] <- na.omit(vi_error_select)
    }

    dd <- as.data.frame(do.call("rbind", vi))
    dd$pixel <- k
    dat[[k]] <- dd
  }

  dat <- do.call("rbind", dat)

  parameters <- data.frame(vi_0 = vi_0_var,
                           vi_delta = vi_delta_var,
                           a = a_var,
                           b = b_var,
                           c = c_var,
                           d = d_var,
                           l = l_var,
                           n = as.numeric(tapply(dat$y, dat$pixel, length)))

  annual_variation <- data.frame(vi_0 = vi_0_var_y,
                                 vi_delta = vi_delta_var_y,
                                 a = a_var_y,
                                 b = b_var_y,
                                 c = c_var_y,
                                 d = d_var_y,
                                 l = l_var_y)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = doy, y = vi, col = y)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~pixel) +
    ggplot2::xlab("\nDay of year") +
    ggplot2::ylab("Vegetation index\n") +
    ggplot2::scale_colour_gradientn("Year", colors = cols) +
    ggplot2::xlim(1, 365)
  print(p)

  # Return everything

  return(list(data = dat, parameters = parameters, annual = annual_variation, plot = p))
}

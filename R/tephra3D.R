Tephra3d = function(vx0, vy0, vz0, x0 = 0, y0 = 0, z0 = 0, t0 = 0, rho_r = 2000, r = 1, dt = 0.01, Cd = 0.6, verbose = FALSE, rho_a = NULL, zt = NULL, wx = 0, wy = 0, wz = 0, mindist = 0, TOPO = NULL, Kh = 0, Kz = 0, eddy_timescale = 60, g = 9.80665){
  #' Tephra Lagrangian Transport Model
  #'
  #' This function models the transport of a single particle through a spatially variable, windy, turbulent atmosphere with gravity. It allows 4D-varying atmospheric wind and density. Currently, only Rayleigh drag (low Re) is permitted.
  #'
  #' @param vx0 initial x component of velocity [m/s]
  #' @param vy0 initial y component of velocity [m/s]
  #' @param vz0 initial z component of velocity [m/s]
  #' @param x0 initial easting [m]
  #' @param y0 initial northing [m]
  #' @param z0 initial elevation [m]
  #' @param t0 initial time [s]
  #' @param rho_r density of tephra particle [kg/m^3]
  #' @param r rock radius [m]
  #' @param dt model time interval [s]
  #' @param Cd drag coefficient [unitless]
  #' @param verbose logical: print runtime info to screen?
  #' @param rho_a scalar or function(x,y,z,t) or function(z) giving atmospheric density [kg/m^3]. If NULL, use a variable-density isothermal atmosphere (T = 0 deg C) 
  #' @param zt function(x,y) giving topographic height [m]
  #' @param wx scalar or function(x,y,z,t) or function(z) giving component of wind to the east [m/s]
  #' @param wy scalar or function(x,y,z,t) or function(z) giving component of wind to the north [m/s]
  #' @param wz scalar or function(x,y,z,t) or function(z) giving upward component of wind [m/s]
  #' @param mindist minimum distance a particle must travel before simulation can stop. This is to prevent early model ends due to spurious collision with crater.
  #' @param TOPO DEM: list containing vectors x and y and matrix z with dimensions (length(x), length(y))
  #' @param Kh horizontal eddy diffusivity [m^2/s]
  #' @param Kz vertical eddy diffusivity (often neglected) [m^2/s]
  #' @param eddy_timescale 1/e decay time of turbulent eddies [s]
  #' @param g gravitational acceleration [m/s^2]
  #' @keywords misc
  #' @export
  #' @examples
  #' Tephra3d(vx0 = 40, vy0 = 0, vz0 = 40, z0 = 0)
  
  ## rock drag coefficient: 0.6 is somewhere between a sphere and low cube

  ## default atmospheric density function: assume sea level pressure at z=0
  ## with fixed temperature profile (273 K) and 6.8 km scale height
  if(is.null(rho_a)){
    T = 273.15 # isothermal atm profile: 0 deg C
    R = 287.04 # gas constant for air
    scale_height = R * T / g 
    rho_a = function(x,y,z,t)101325/(R * T) * exp(-z/scale_height) 
  }

  ## default topography function: assumes flat topo with z = 0
  if(is.null(zt)){
    zt = function(x,y) x*0
  }

  ## default topography function if TOPO DEM is provided
  if(!is.null(TOPO)){
    zt = function(x,y)TopoInterp(x, y, TOPO)
  }

  ## check all the atmospheric functions to see how many arguments they accept.
  ## if just one, assume that they vary with z only, and add a wrapper for compatibility.
  ## if one of the functions is actually just a scalar, make it into a function
  if(is.function(rho_a) && length(formals(rho_a)) == 1){ # atm density function
    rho_a_inner = rho_a
    rho_a = function(x,y,z,t)rho_a_inner(z)
  }else if(is.numeric(rho_a)){
    rho_a_inner = rho_a
    rho_a = function(x,y,z,t)rho_a_inner
  }    
  if(is.function(wx) && length(formals(wx)) == 1){ # zonal wind function
    wx_inner = wx
    wx = function(x,y,z,t)wx_inner(z)
  }else if(is.numeric(wx)){
    wx_inner = wx
    wx = function(x,y,z,t)wx_inner
  }    
  if(is.function(wy) && length(formals(wy)) == 1){ # meridional wind function
    wy_inner = wy
    wy = function(x,y,z,t)wy_inner(z)
  }else if(is.numeric(wy)){
    wy_inner = wy
    wy = function(x,y,z,t)wy_inner
  }    
  if(is.function(wz) && length(formals(wz)) == 1){ # vertical wind function
    wz_inner = wz
    wz = function(x,y,z,t)wz_inner(z)
  }else if(is.numeric(wz)){
    wz_inner = wz
    wz = function(x,y,z,t)wz_inner
  }    

  m = 4/3*pi*r^3 * rho_r # tephra particle mass
  A = pi*r^2 # cross-sectional area of sphere

  ## initialize position, velocity, time vectors
  vx = vx0
  vy = vy0
  vz = vz0
  x = x0
  y = y0
  z = z0
  t = t0

  ## initialize turbulent air velocity
  tx = 0
  ty = 0
  tz = 0

  ## loop over time until the particle hits the ground
  i = 1
  while(z[i] >= zt(x[i],y[i]) || (abs(x[i]) < mindist)){
    i = i+1
    ## these following equations taken from Maeno et al., 2013: Ballistic ejecta and eruption condition of the vulcanian explosion of Shinmoedake volcano, Kyushu, Japan on 1 February, 2011, and modified for wind only.  The v_x * v is necessary in order to preserve magnitude of drag force independent of coordinate system, even though it looks weird that f_x depends on v_z.  This can be explained by the mass of air having to be displaced by the projectile increases with the velocity magnitude, and that mass affects both components of the drag force.

    ## determine total air velocity: mean wind plus turbulent
    ##   turbulent air velocity: v[i+1] = v[i] * (1 - dt/tau) + sqrt(2 K dt)/tau * r(0,1)
    tx = tx * (1 - dt/eddy_timescale) + sqrt(2*Kh*dt)/eddy_timescale * rnorm(1, 0, 1)
    ty = ty * (1 - dt/eddy_timescale) + sqrt(2*Kh*dt)/eddy_timescale * rnorm(1, 0, 1)
    tz = tz * (1 - dt/eddy_timescale) + sqrt(2*Kz*dt)/eddy_timescale * rnorm(1, 0, 1)

    ## wtx, wty, wtz are total air velocities (mean + eddy)
    wtx = wx(x = x[i-1], y = y[i-1], z = z[i-1], t = t[i-1]) + tx
    wty = wy(x = x[i-1], y = y[i-1], z = z[i-1], t = t[i-1]) + ty
    wtz = wz(x = x[i-1], y = y[i-1], z = z[i-1], t = t[i-1]) + tz

    ## calculate new particle velocities
    v_mag = sqrt((wtx - vx[i-1])^2+(wty - vy[i-1])^2 + (wtz-vz[i-1])^2) + 1e-12 # have to add an epsilon to avoid divide-by-zero errors later

    ## calculate local air density
    rho_i= rho_a(x = x[i-1], y = y[i-1], z = z[i-1], t = t[i-1])
    ## magnitude of drag acceleration, to be applied proportionally to velocity components 
    drag = 0.5 * rho_i* v_mag^2 * A * Cd / m

    ## apply drag force to horizontal components
    vx[i] = vx[i-1] - (vx[i-1] - wtx)/v_mag * drag * dt
    vy[i] = vy[i-1] - (vy[i-1] - wty)/v_mag * drag * dt

    ## apply drag force and gravity-buoyancy force to vertical component
    vz[i] = vz[i-1] - ((vz[i-1] - wtz)/v_mag * drag) * dt - (g * (rho_r - rho_i)/rho_r) * dt
    z[i] = z[i-1] + (vz[i] + vz[i-1])/2 * dt
    x[i] = x[i-1] + (vx[i] + vx[i-1])/2 * dt
    y[i] = y[i-1] + (vy[i] + vy[i-1])/2 * dt
    t[i] = t[i-1] + dt
    if(verbose){
      print(paste(i, z[i], x[i], y[i]))
    }
    if(is.na(zt(x[i],y[i])) ){
      warning('Out of bounds: topographic surface undefined for x or y')
      z[i] = NaN
      break
    }
  }

  ## if the loop breaks on the first iteration, it means the particle starts underground. this should be a meaningless NaN result.
  if(i == 1){
    x = c(x,NaN)
    y = c(y,NaN)
    z = c(z,NaN)
    t = c(t,NaN)
  }
  ## After leaving the loop, the particle is now below the interpolated topography.
  ## Now, interpolate to find where the particle would actually hit the topography.
  l = sqrt( (x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 )
  slope_tephra = (z[i] - z[i-1]) / l # slope of tephra path
  slope_topo = ( zt(x[i], y[i]) - zt(x[i-1], y[i-1]) ) / l # slope of topo
  l_intersect = (z[i-1] - zt(x[i-1], y[i-1])) / (slope_topo - slope_tephra) # point along (x[i-1], y[i-1])-(x[i], y[i]) line where topo and tephra meet--linear approximation

  ## now, replace x[i], y[i], z[i], and t[i] with interpolated values.
  ## Note that this can have small errors if either topography or trajectories are significantly curved.
  x[i] = x[i-1] + (x[i] - x[i-1]) * l_intersect / l
  y[i] = y[i-1] + (y[i] - y[i-1]) * l_intersect / l
  t[i] = t[i-1] + (t[i] - t[i-1]) * l_intersect / l 
  z[i] = zt(x[i], y[i]) ## force z[i] to be on the Gaussian-interpolated topographic surface
  
  return(list(x = x, y = y, z = z, t = t, xf = x[i], yf = y[i], zf = z[i], tf = t[i], vx = vx, vy = vy, vz = vz, vxf = vx[i], vyf = vy[i], vzf = vz[i]))
}



########################################
########################################
########################################


BlastSim3d = function(v, th_i = 2* 1:40, th_a = 0, dt = 0.01, ...){
  #' Explosive Tephra Dispersion Model
  #'
  #' Models the transport of particles ejected at the same velocity and different angles using tephra3d.
  #'
  #' @param v initial velocity [m/s]
  #' @param th_i initial inclination angles to test [degrees]
  #' @param th_a initial azimuth angles to test [degrees clockwise from north]
  #' @param dt model time interval [s]
  #' @param ... parameters to pass to Tephra3d
  #' @keywords misc
  #' @export
  #' @examples
  #' BlastSim3d(v = 10, th_i = 2* 1:40, th_a = 0, dt = 0.01)
    cosd = function(x)cos(pi/180*x)
    sind = function(x)sin(pi/180*x)
    L = meshgridn(list(th_i, th_a)); names(L) = c('th_i', 'th_a')
    L = c(L, v = v, dt = dt, x = list(), z = list())
    for(i in 1:length(L$th_i)){
        print(c(i, sind(L$th_i[i])*cosd(L$th_a[i])*v))
        T = Tephra3d(vx0 = sind(L$th_i[i])*sind(L$th_a[i])*v, vy0 = sind(L$th_i[i])*cosd(L$th_a[i])*v, vz0 = cosd(L$th_i[i]) * v, dt = dt, ...)
        L$x[[i]] = T$x
        L$y[[i]] = T$y
        L$z[[i]] = T$z
    }
    invisible(L)
}


###########


BlastAnim3d = function(L, tframe = 0.1, dir = '.', TOPO = NULL, az = 0, xlim = NULL, ylim = NULL, zlim = NULL, units = 'm', plotMapView = TRUE, plotCrossSection = TRUE){
  #' Tephra Transport Snapshots
  #'
  #' Generates png files showing map view and cross-section view of tephra motion and final position on ground
  #'
  #' @param L output of BlastSim3d
  #' @param tframe time interval between frames [s]
  #' @param dir directory where png files should be saved
  #' @param TOPO DEM: list containing vectors x and y and matrix z with dimensions (length(x), length(y))
  #' @param az azimuth of section line (degrees clockwise from North)
  #' @param xlim easting limits for map view 
  #' @param ylim northing limits for map view
  #' @param zlim elevation limits for section view
  #' @param units units of length (string)
  #' @param plotMapView logical: should the map view panel be plotted?
  #' @param plotCrossSection logical: should the cross-section panel be plotted?
  #' @keywords misc
  #' @export
  #' @import graphics
  #' @import grDevices
  #' @import stats
  #' @examples
  #' ## Not run:
  #' ## BlastSim3d(v = 40, th_i = 2* 1:40, th_a = 90, dt = 0.01)
  #' ## BlastAnim3d(L, tframe = 0.1, az = 90)
  #' ## ImageMagick shell command: animate -delay 10 * # animate with 0.1-s frame rate
  tephraColors = colorRampPalette(c('turquoise', 'darkblue', 'purple', 'red', 'orange'), bias = 1, space = c("rgb"))
  ## define zt, function to interpolate digital elevation model
  if(!is.null(TOPO)){
    zt = function(l)sapply(l, function(x)TopoInterp(sin(pi/180*az)*x, cos(pi/180*az)*x, TOPO))
  }else{
    zt = function(l)0*l
  }
  
  if(is.null(xlim)){
    xlim = range(unlist(L$x), na.rm = TRUE)
    xlim = xlim + c(-1, 1) * 0.25 * diff(xlim)
  }
  if(is.null(ylim)){
    ylim = range(unlist(L$y), na.rm = TRUE)
    ylim = ylim + c(-1, 1) * 0.25 * diff(ylim)
  }
  if(is.null(zlim)){
    zlim = range(unlist(L$z), na.rm = TRUE)
    zlim = zlim + c(-1, 1) * 0.25 * diff(zlim)
  }
  llim = sin(az*pi/180)*xlim + cos(az*pi/180)*ylim
  rK = function(K, L)if(!is.null(K))(K/max(L$K))^0.25 else 1
  ##    rz = function(z0, L)if(!is.null(z0))rgb(1-(z0 - min(L$z0))/diff(range(L$z0)), 0, (z0 - min(L$z0))/diff(range(L$z0))) else 'red'
  ##rz = function(K, L)if(!is.null(K))rgb((log(K) - min(log(L$K)))/diff(range(log(L$K))), 0, 1-(log(K) - min(log(L$K)))/diff(range(log(L$K)))) else 'red'
  rz = function(K, L)if(!is.null(K)) tephraColors(100)[ceiling(100 * (log(K) - min(log(L$K)))/diff(range(log(L$K))))] else 'red'
                                        #n = ceiling((5 + max(sapply(L$x, length)) * L$dt) / tframe) # add 5 s after the last ballistic lands
  n = ceiling(max(sapply(L$x, length)) * L$dt / tframe) 

  nc = 1 + floor(log(n, 10))
  if(is.null(L$K)){
    w = 1:length(L$x)
  }else{
    w = order(L$K, decreasing = TRUE)
  }
  for(i in 1:n){
    print(paste(i, 'of', n))
    k = i * tframe/L$dt
    nn = substr(10^nc + i, 2, nc + 1)
    png(paste(dir, '/', nn, '.png', sep = ''), width = 960, height = 480, pointsize = 18)
###### l-z plot
    if(((i * tframe) %% 1) == 0){
      main = paste('t = ', i*tframe, '.0 s', sep = '')
    }else{
      main = paste('t =', i*tframe, 's')
    }
    par(mgp = c(1.75, 0.5, 0), mar = c(2.75, 2.75, 2.75, 1.5), mfrow = c(1, plotMapView + plotCrossSection))
    if(plotCrossSection){
      plot(NaN, xlim = llim, ylim = zlim, main = main, xaxt='n', yaxt='n', ylab = '', xlab = '')
      if(units == 'm'){
        title(xlab = 'Radial Distance (m)', ylab = 'Elevation (m ASL)')
        axis(side = 1, at = pretty(llim), labels = pretty(llim))
        axis(side = 2, at = pretty(zlim), labels = pretty(zlim))
      }else if(units == 'km'){
        title(xlab = 'Radial Distance (km)', ylab = 'Elevation (km ASL)')
        axis(side = 1, at = pretty(llim), labels = pretty(llim)/1000)
        axis(side = 2, at = pretty(zlim), labels = pretty(zlim)/1000)
      }            
      for(j in w){
        if(is.na(L$x[[j]][k])){
          wm = max(which(!is.na(L$x[[j]])))
          points(sin(pi/180*az)*L$x[[j]][wm] + cos(pi/180*az)*L$y[[j]][wm], L$z[[j]][wm], col = rz(L$K[j], L), pch = 19, cex = rK(L$K[j], L))
        }else{
          points(sin(pi/180*az)*L$x[[j]][k] + cos(pi/180*az)*L$y[[j]][k], L$z[[j]][k], col = rz(L$K[j], L), pch = 19, cex = rK(L$K[j], L))
        }
      }
      if(!is.null(zt)){
        lines(seq(llim[1], llim[2], length.out = 500), zt(seq(llim[1], llim[2], length.out = 500)))
      }
    }# if plotCrossSection
############ x-y plot
    if(plotMapView){
      plot(NaN, xlim = xlim, ylim = ylim, asp = 1, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
      if(!is.null(TOPO)){
        image(TOPO$x, TOPO$y, TOPO$z, col = terrain.colors(30), add = TRUE)
        if(units == 'm'){
          title(xlab = 'East (m)', ylab = 'North (m)')
          axis(side = 1, at = pretty(xlim), labels = pretty(xlim))
          axis(side = 2, at = pretty(ylim), labels = pretty(ylim))
        }else if(units == 'km'){
          title(xlab = 'East (km)', ylab = 'North (km)')
          axis(side = 1, at = pretty(xlim), labels = pretty(xlim)/1000)
          axis(side = 2, at = pretty(ylim), labels = pretty(ylim)/1000)
        }            
        wx = which(TOPO$x >= xlim[1] & TOPO$x <= xlim[2])
        wy = which(TOPO$y >= ylim[1] & TOPO$y <= ylim[2])
        levels = pretty(range(TOPO$z[wx,wy]), 10)
        contour(TOPO$x, TOPO$y, TOPO$z, col = terrain.colors(30), drawlabels=FALSE, add = TRUE, levels = levels)
        abline(a=0, b=1/tan(az*pi/180+1e-6)) # plot the section line: +1e-6 so that azimuth zero still works.
      }
      for(j in w){
        if(is.na(L$x[[j]][k])){
          wm = max(which(!is.na(L$x[[j]])))
          points(L$x[[j]][wm], L$y[[j]][wm], col = rz(L$K[j], L), pch = 19, cex = rK(L$K[j], L))
        }else{
          points(L$x[[j]][k], L$y[[j]][k], col = rz(L$K[j], L), pch = 19, cex = rK(L$K[j], L))
        }
      }
    } # if plotMapView
    dev.off()
  }
}

bracket = function(x, xlist){
    which((x - xlist) %in% c(max(x - xlist[xlist > x]), min(x - xlist[xlist < x])))
}
TopoInterp = function(x, y, TOPO = NULL, N = 10){
  #' Topography Interpolation
  #'
  #' Interpolates elevation at point (x, y) given Digital Elevation Model (DEM).
  #'
  #' @param x Easting of point to interpolate [same units as TOPO$x, TOPO$y]
  #' @param y Northing of point to interpolate [same units as TOPO$x, TOPO$y]
  #' @param TOPO DEM: list containing vectors x and y and matrix z with dimensions (length(x), length(y))
  #' @param N Smoothing parameter, must be positive. Larger N means less smoothing. 
  #' @keywords misc
  #' @export
  #' @examples
  #' data(VILL)
  #' contour(VILL, xlim = c(-500, 500), ylim = c(-500, 500))
  #' TopoInterp(0, 0, VILL) # interpolate elevation at point (0, 0)
 
    # scalar input only
    wx = bracket(x, TOPO$x)
    wy = bracket(y, TOPO$y)
    # super ugly code to make the window 6 long instead of 2
    if(min(wx) != 1) wx = c(min(wx)-1, wx)
    if(min(wx) != 1) wx = c(min(wx)-1, wx)
    if(max(wx) != length(TOPO$x)) wx = c(wx, max(wx)+1)
    if(max(wx) != length(TOPO$x)) wx = c(wx, max(wx)+1)
    if(min(wy) != 1) wy = c(min(wy)-1, wy)
    if(min(wy) != 1) wy = c(min(wy)-1, wy)
    if(max(wy) != length(TOPO$y)) wy = c(wy, max(wy)+1)
    if(max(wy) != length(TOPO$y)) wy = c(wy, max(wy)+1)
    
    dx = abs(x - TOPO$x[wx])
    dy = abs(y - TOPO$y[wy])

    R = sqrt(matrix(dx, length(dx), length(dy))^2 + matrix(dy, length(dx), length(dy), byrow = TRUE)^2)
    M = exp(-N*R^2/mean(R^2))/sum(exp(-N*R^2/mean(R^2))) # gaussian--N = 10 works well

    sum(M * TOPO$z[wx,wy])
}

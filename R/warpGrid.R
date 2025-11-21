#     warpGrid.R Warp grid to allow plotting in a different projection
#
#     Copyright (C) 2019 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Grid warping
#' @description Warp grid to allow plotting in a different projection.
#' @param data A C4R grid (or multimember C4R grid) object, or climatology C4R grid.
#' @param original.CRS character or object as passed to function \code{\link[sf]{st_crs}} with the original projection.
#' Default to longlat projection (\code{"+init=epsg:4326"}).
#' @param new.CRS character string or object, as passed to function \code{\link[sf]{st_crs}}, specifying the target projection.
#' Default to polar stereographic projection (\code{"+init=epsg:3995"}).
#' @param int.method Resampling method. Default to \code{"bilinear"}. See details.
#' 
#' @return Warped grid with the structure of a C4R grid.
#' 
#' @details 
#' This function is a wrapper of the gdal warping capabilities via \code{stars::st_warp}.
#' 
#'  \strong{int.method}
#'  
#'  By default bilinear interpolation is applied to get a complete grid in the target projection. Other options are \code{"near"}, \code{"cubic"},
#'   \code{"cubicspline"} etc., passed to the argument \code{method} in \code{stars::st_warp}.

#' @export
#' @importFrom sf st_crs
#' @importFrom stars st_as_stars st_warp
#' @import transformeR
#' @author A. Casanueva, J. Bedia, M. Iturbide
#' @examples
#' library(climate4R.datasets)
#' data(ncep_hgt500_2000)
#' grid <- warpGrid(climatology(ncep_hgt500_2000))
#' # Example of application: plot in polar stereographic projection
#' library(visualizeR)
#' l1 <- get(load(paste0(find.package("visualizeR"), "/countries.rda"))) # world coastline
#' l1 <- sf::st_transform(sf::st_as_sf(l1[[2]]), crs = attr(grid$xyCoords, "projection"))
#' visualizeR::spatialPlot(grid, sp.layout = list(list(l1, first = FALSE)))

warpGrid <- function(data,
                     original.CRS = "+init=epsg:4326",
                     new.CRS = "+init=epsg:3995", 
                     int.method = "bilinear") {
  
  # *** Convert grid to sp ***
  pattern <- transformeR::grid2sp(data) 

  # *** Convert sp to stars ***
  suppressWarnings(sp::proj4string(pattern) <- sp::CRS(NA_character_)) # Remove invalid CRS from sp so that stars can handle it
  pattern_stars <- stars::st_as_stars(pattern) # Convert to stars
  sf::st_crs(pattern_stars) <- sf::st_crs(original.CRS) # Assign valid CRS to the stars object
  band_names <- names(pattern_stars) # Preserve members 

  # *** IMAGE RE-PROJECTION ***
  warped_list <- lapply(band_names, function(nm) {
    stars::st_warp(
      pattern_stars[nm], 
      dest     = pattern_stars[nm],
      crs      = sf::st_crs(new.CRS),  
      method   = int.method,
      use_gdal = TRUE
    )
  })

  # *** Convert stars to sp ***

  # First member
  warped_sp <- as(warped_list[[1]], "Spatial")

  # Add extra members as columns 
  if (length(warped_list) > 1) { 
    extra_cols <- lapply(warped_list[-1], function(w) {
      sp_tmp <- as(w, "Spatial")
      sp_tmp@data[, 1]         
    })
    warped_sp@data <- data.frame(
      warped_sp@data,
      do.call(cbind, extra_cols)
    )
  }

  # Make sure column names match original member names
  names(warped_sp@data) <- band_names

  outf <- newf <- NULL
  
  # *** Convert sp to grid ***
  start <- getRefDates(data, which = "start")
  end <- getRefDates(data, which = "end")
  
  grid <- transformeR::sgdf2clim(sp = warped_sp,
                                 varName = getVarNames(data),
                                 level = getGridVerticalLevels(data),
                                 dates = list(start = start, end = end),
                                 season = getSeason(data))
}


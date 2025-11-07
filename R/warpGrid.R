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
#' This function is a wrapper of the gdal warping capabilities via gdalwarp.  
#' 
#'  \strong{int.method}
#'  
#'  By default bilinear interpolation is applied to get a complete grid in the target projection. Other options are \code{"near"}, \code{"cubic"},
#'   \code{"cubicspline"} etc., passed to the argument \code{r} in \code{gdalUtils::gdalwarp}.

#' @export
#' @importFrom sf st_as_sf st_coordinates st_crs st_drop_geometry
#' @importFrom gdalUtils gdalwarp
#' @importFrom stars st_as_stars read_stars write_stars
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
  
  # *** Check for members ***
  nmem <- getShape(data, "member")
  member <- ifelse((nmem == 1 | is.na(nmem)), FALSE, TRUE)
  
  crs_from_input <- function(value, arg_name) {
    if (inherits(value, "crs")) {
      return(value)
    }
    if (is.null(value) || (length(value) == 1 && is.na(value))) {
      return(sf::st_crs(NA))
    }
    tryCatch({
      sf::st_crs(value)
    }, error = function(err) {
      stop("Non-valid ", arg_name, " argument")
    })
  }
  crs_to_string <- function(crs_obj) {
    if (is.na(crs_obj)) {
      return(NA_character_)
    }
    if (!is.null(crs_obj$input) && !is.na(crs_obj$input)) {
      return(crs_obj$input)
    }
    if (!is.null(crs_obj$wkt) && !is.na(crs_obj$wkt)) {
      return(crs_obj$wkt)
    }
    NA_character_
  }
  original_crs <- crs_from_input(original.CRS, "original.CRS")
  new_crs <- crs_from_input(new.CRS, "new.CRS")
  original_crs_string <- crs_to_string(original_crs)
  new_crs_string <- crs_to_string(new_crs)

  # *** CONVERT GRID TO STARS ***
  pattern <- transformeR::grid2sp(data)
  pattern_stars <- stars::st_as_stars(pattern)
  sf::st_crs(pattern_stars) <- original_crs

  # *** WRITE A GDAL GRID MAP ***
  outf <- tempfile(fileext = ".tif")
  suppressWarnings(
    stars::write_stars(pattern_stars, dsn = outf, driver = "GTiff", NA_value = NA_real_)
  )

  # *** IMAGE RE-PROJECTION ***
  newf <- tempfile(fileext = ".tif")
  s_srs_value <- if (is.na(original_crs_string)) NULL else original_crs_string
  t_srs_value <- if (is.na(new_crs_string)) NULL else new_crs_string
  suppressMessages(
    gdalUtils::gdalwarp(srcfile = outf,
                        s_srs = s_srs_value,
                        t_srs = t_srs_value,
                        dstfile = newf,
                        r = int.method)
  )
  # *** READ NEW IMAGE ***
  warped <- stars::read_stars(newf, NA_value = NA_real_)
  outf <- newf <- NULL

  # *** stars2grid ***
  start <- getRefDates(data, which = "start")
  end <- getRefDates(data, which = "end")

  warped_sf <- sf::st_as_sf(warped, as_points = TRUE, merge = TRUE)
  coords_matrix <- sf::st_coordinates(warped_sf)
  coords_df <- data.frame(x = coords_matrix[, 1], y = coords_matrix[, 2])
  order_idx <- order(coords_df$y, coords_df$x)
  coords_df <- coords_df[order_idx, , drop = FALSE]
  values_matrix <- as.matrix(sf::st_drop_geometry(warped_sf))[order_idx, , drop = FALSE]
  values_matrix[is.nan(values_matrix)] <- NA
  data_matrix <- t(values_matrix)
  x_vals <- unique(coords_df$x)
  y_vals <- unique(coords_df$y)
  gridded_data <- mat2Dto3Darray(data_matrix, x_vals, y_vals)
  attr(gridded_data, "dimensions") <- attr(data$Data, "dimensions")

  grid <- data
  grid$Data <- gridded_data
  if (is.list(grid$xyCoords)) {
    grid$xyCoords$x <- x_vals
    grid$xyCoords$y <- y_vals
    attr(grid$xyCoords, "resX") <- if (length(x_vals) > 1) x_vals[2] - x_vals[1] else 0
    attr(grid$xyCoords, "resY") <- if (length(y_vals) > 1) y_vals[2] - y_vals[1] else 0
  } else {
    grid$xyCoords <- coords_df
    attr(grid$xyCoords, "resX") <- 0
    attr(grid$xyCoords, "resY") <- 0
  }
  attr(grid$xyCoords, "projection") <- new_crs_string

  attr(grid$Data, "dimensions") <- attr(data$Data, "dimensions")
  grid$Dates <- list(start = start, end = end)
  grid$Season <- getSeason(data)
  return(grid)
}


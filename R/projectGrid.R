#     projectGrid.R Grid datum definition and transformation
#
#     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
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


#' @title grid datum definition and transformation
#' @description Defines and/or transforms the projection of a grid (or station data) by means of an \code{\link[sf]{st_crs}} object.
#' @param grid a grid or multigrid (including station data).
#' @param original.CRS character or object as passed to function \code{\link[sf]{st_crs}}. It defines the original projection. If the data contains the
#' projection information, a warning is returned and the projection in redefined.
#' @param new.CRS character or object as passed to function \code{\link[sf]{st_crs}}.
#' @details This function uses \code{\link[sf]{st_transform}} and \code{\link[sf]{st_crs}} from package \pkg{sf}.
#' @seealso \code{\link[sf]{st_transform}}, \code{\link[sf]{st_crs}}.
#' 
#' @author M. Iturbide
#' @export
#' @importFrom sf st_as_sf st_coordinates st_crs st_transform
#' @importFrom abind abind
#' @import transformeR
#' @examples
#' library(climate4R.datasets)
#' data("VALUE_Iberia_pr")
#' plot(getCoordinates(VALUE_Iberia_pr))
#' grid <- projectGrid(VALUE_Iberia_pr,
#'                     original.CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
#'                     new.CRS = "+init=epsg:28992")
#' plot(getCoordinates(grid))
#' 
#' data("EOBS_Iberia_pr")
#' plot(get2DmatCoordinates(EOBS_Iberia_pr))
#' grid <- projectGrid(EOBS_Iberia_pr,
#'                     original.CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
#'                     new.CRS = "+init=epsg:28992")
#' plot(get2DmatCoordinates(grid))
#' require(visualizeR)
#' spatialPlot(climatology(grid))




projectGrid <- function(grid,
                        original.CRS = "",
                        new.CRS = "") {
  orig.datum <- attr(grid$xyCoords, "projection")
  # if (orig.datum == "RotatedPole") stop("This function is not applicable to this projection. See Details")
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
  if (!is.null(orig.datum) & !is.na(original_crs)) {
    warning("CAUTION! Grid with previusly defined projection: ", orig.datum)
    attr(grid$xyCoords, "projection") <- crs_to_string(original_crs)
  } else if (is.null(orig.datum) & !is.na(original_crs)) {
    attr(grid$xyCoords, "projection") <- crs_to_string(original_crs)
  } else if (is.null(orig.datum) & is.na(original_crs)) {
    stop("Please define original.CRS")
  } else if (!is.null(orig.datum) & is.na(original_crs)) {
    original_crs <- tryCatch({
      crs_from_input(orig.datum, "original.CRS")
    }, error = function(err) {
      stop("Grid with non-valid defined projection. Please, use argument original.CRS to redefine it correctly")
    })
    attr(grid$xyCoords, "projection") <- crs_to_string(original_crs)
  }
  if (is.null(attr(grid$xyCoords, "projection"))) {
    attr(grid$xyCoords, "projection") <- crs_to_string(original_crs)
  }
  data <- get2DmatCoordinates(grid)
  coords_df <- as.data.frame(data)
  if (ncol(coords_df) >= 2) {
    colnames(coords_df)[1:2] <- c("x", "y")
  }
  sppoints <- sf::st_as_sf(coords_df, coords = c("x", "y"), crs = original_crs)
  message("[", Sys.time(), "] ", "Arguments of the original projection defined as ", crs_to_string(original_crs))
  if (!is.na(new_crs)) {
    message("[", Sys.time(), "] ", "Projecting..")
    sppoints.new <- sf::st_transform(sppoints, new_crs)
    new.coords <- sf::st_coordinates(sppoints.new)
    x <- unique(new.coords[,1])
    y <- unique(new.coords[,2])
    if (length(x) > 1 & length(y) > 1) { # a single location?
      xdists <- lapply(1:(length(x) - 1), function(l) {
        x[l + 1] - x[l]
      })
      ydists <- lapply(1:(length(y) - 1), function(l) {
        y[l + 1] - y[l]
      })
      xa <- sum(unlist(xdists) - unlist(xdists)[1])
      ya <- sum(unlist(ydists) - unlist(ydists)[1])
      cond <- any(abs(c(xa, ya)) > 1e-05) # regular coordinates?
    } else {
      cond <- TRUE # a single location considered as irregular
    }
    if (cond) {
      if (isRegular(grid)) { 
        grid <- redim(grid, member = TRUE, runtime = TRUE)
        data.aux1 <- lapply(1:getShape(grid)["runtime"], function(r) {
          data.aux0 <- lapply(1:getShape(grid)["member"], function(m) {
            array3Dto2Dmat(redim(subsetGrid(grid, runtime = r, members = m), member = FALSE)$Data)
          })
          do.call("abind", list(data.aux0, along = 0))
        })
        grid$Data <-  do.call("abind", list(data.aux1, along = 0))
        attr(grid$Data, "dimensions") <- c("runtime", "member", "time", "loc")
        grid <- redim(redim(grid, drop = T), member = FALSE, loc = TRUE)
      }
      grid$xyCoords <- as.data.frame(new.coords)
      colnames(grid$xyCoords) <- c("x", "y")
      attr(grid$xyCoords, "projection") <- crs_to_string(new_crs)
      attr(grid$xyCoords, "resX") <- 0
      attr(grid$xyCoords, "resY") <- 0
    } else {
      grid$xyCoords <- list("x" = unique(new.coords[,1]), "y" = unique(new.coords[,2]))
      attr(grid$xyCoords, "projection") <- crs_to_string(new_crs)
      attr(grid$xyCoords, "resX") <- xdists[[1]]
      attr(grid$xyCoords, "resY") <- ydists[[1]]
    }
    message("[",Sys.time(), "] ", "Done.")
  }
  return(grid)
}

#end

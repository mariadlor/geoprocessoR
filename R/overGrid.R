##     overGrid.R Grid scaling
##
##     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program. If not, see <http://www.gnu.org/licenses/>.


#' @title Spatial overlay for grids and spatial objects
#' @description Application of spatial overlay using \pkg{sf} for grids or
#' station data.
#' @param grid Input grid or station data.
#' @param layer Spatial object from which the geometries or attributes are
#' queried (see \code{help(package = "sf")}.
#' @param subset Logical (default is FALSE). If TRUE, spatial subsetting is performed over
#' the imput grid. Otherwise, grid-boxes outside \code{layer} are assigned as NA.
#' 
#' @details All grid locations outside layer are filled with NAs. 
#' 
#' @return A grid or station data.
#'
#' @importFrom sf st_as_sf st_bbox st_coordinates st_crs st_drop_geometry st_intersects
#' @import transformeR
#' 
#' @return A grid
#' @author M. Iturbide
#' @family subsetting
#' @export

overGrid <- function(grid, layer, subset = FALSE) {
      # commented lines are for posible runtime dimension consideration,
      # in which case function bindGrid.runtime should be created
      # grid <- redim(grid, runtime = TRUE)
      layer_sf <- if (inherits(layer, "sf")) layer else sf::st_as_sf(layer)
      if (subset) {
            layer_bbox <- sf::st_bbox(layer_sf)
            grid <- subsetGrid(grid,
                               lonLim = c(layer_bbox["xmin"], layer_bbox["xmax"]),
                               latLim = c(layer_bbox["ymin"], layer_bbox["ymax"]),
                               outside = TRUE)
      }
      loc <- "loc" %in% getDim(grid)
      grid <- redim(grid, loc = loc)
      grid_crs_attr <- attr(grid$xyCoords, "projection")
      grid_crs <- tryCatch({
            if (is.null(grid_crs_attr) || is.na(grid_crs_attr)) {
                  sf::st_crs(NA)
            } else {
                  sf::st_crs(grid_crs_attr)
            }
      }, error = function(err) sf::st_crs(NA))
      layer_crs <- sf::st_crs(layer_sf)
      if (!is.na(layer_crs) && !is.na(grid_crs) && !identical(layer_crs$wkt, grid_crs$wkt)) {
            layer_sf <- sf::st_transform(layer_sf, grid_crs)
      }
      # n.run <- getShape(grid)["runtime"]
      n.mem <- getShape(grid)["member"]
      if (loc) {
            coords <- grid$xyCoords[, 2:1]
      } else {
            coords <- expand.grid(getCoordinates(grid)$y, getCoordinates(grid)$x)
      }
      # grr <- lapply(1:n.run, function(k){
            # grid.r <- subsetGrid(grid, runtime = k)
            grm <- lapply(1:n.mem, function(x) {
                  grl <- redim(subsetGrid(grid, members = x), member = FALSE, loc = loc)
                  dimNames.sub <- getDim(grl)
                  if (loc) {
                        dat <- grl$Data
                  } else {
                        dat <- array3Dto2Dmat(grl$Data)
                  }
                  values_df <- data.frame(t(dat))
                  values_df$..x <- coords[, 2]
                  values_df$..y <- coords[, 1]
                  points_sf <- sf::st_as_sf(values_df, coords = c("..x", "..y"), crs = grid_crs)
                  inside <- lengths(sf::st_intersects(points_sf, layer_sf)) > 0
                  data_cols <- setdiff(names(points_sf), attr(points_sf, "sf_column"))
                  if (length(data_cols) > 0) {
                        points_sf[!inside, data_cols] <- NA
                  }
                  if (subset & loc) {
                        points_sf <- points_sf[inside, , drop = FALSE]
                  }
                  coords_sf <- as.data.frame(sf::st_coordinates(points_sf))
                  colnames(coords_sf) <- c("x", "y")
                  data_sf <- sf::st_drop_geometry(points_sf)
                  if (loc) {
                        grl$Data <- unname(as.matrix(t(data_sf)))
                        grl$xyCoords <- coords_sf
                  } else {
                        grl$Data <- mat2Dto3Darray(t(as.matrix(data_sf)), unique(coords_sf$x), unique(coords_sf$y))
                        grl$xyCoords$x <- unique(coords_sf$x)
                        grl$xyCoords$y <- unique(coords_sf$y)
                  }
                  attr(grl$Data, "dimensions") <- dimNames.sub
                  grl
            })
            if (n.mem > 1) {
                  newgrid <- do.call("bindGrid", c(grm, dimension = "member"))
            } else {
                  newgrid <- grm[[1]]
            }
      # })
        # newgrid <-  do.call("bindGrid.runtime", grr)
      return(newgrid)
}





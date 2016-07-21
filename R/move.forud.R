#' A For Loop creating Utilization Distribution Contours for dbbmm (Dynamic Brownian Bridge Movement Model)
#'
#'\code{move.forud} function that creates step utilization distribution contours based on the \code{brownian
#'.bridge.dyn}from package move and exports the created contour lines as a \code{SpatialPolygonsDataFrame} to a #'set of shapefiles for further whatever in ArcMap (or GIS program of choice).
#'
#'@details  This function is a different approach than \code{move.contour} (which will be deprecated in 2016)
#'that creates the same output (shapefiles of utilization distribution polygons) based on the dbbmm object.
#'The primary difference between \code{move.contour} and \code{move.forud} is that \code{move.forud}
#'uses a \code{for} loop to loop through the range of line segments for which it creates the UD's and then
#'exporting the output as shapefiles one at time to the user-defined path. The time it takes to do each method
#'seems to be the same as the slowness of the process is due to the calculations being done by
#'\code{brownian.bridge.dyn()} of the UD over the raster extent. \code{move.forud} does not use all the memory in the system #'as did \code{move.contour} which saved the objects for one batch export of the shapefiles. \code{move.forud} #'is not particularly any more efficient than \code{move.contour}, but if you have a large raster extent and #'many individual time step utilization distributions that are being estimated, \code{move.forud}
#'does not stack all the raster files in memory, causing issues.
#'
#'I thought that vectorisation would speed this up, but the vast majority of the time is spend in the call to #'\code{brownian.bridge.dyn()} and not the rest of the for loop.  So, I am working on some ways to reduce the
#'background raster on which the contour estimates are created as there are alot of them that are way outside #'the range of useful at the individual time-step level, but when you look at the entire period which \code{brow#'nian.bridge.dyn} does, are actually useful. Think of it as pulling one hour worth of movements out of a year #'worth of GPS locations, it is likely that the 50m in that hour is a small segment of the overall range and #'hence background raster size.
#'
#' As an aside, I do know that it is better to access values in slots using slot(object, 'name') versus the '@@' #'sign, but, the routine I used to pull the requisite values from the DDBMvar object would not work as the data #'necessary does not have slot values assigned to it at this time.
#'
#' @name move.forud
#' @param x, dbbmm object created using \code{brownian.bridge.dyn} function call from \code{move}
#' @param range.subset, range of model step segments that the user wants to create the UD contours for.  Has to be greater than the margin based on the dbbmm object as estimates of variance for
#' the first few steps are not estimable depending on how big of a margin is specified.  Easiest way to define #' this range is to use the row number from the input file.
#' @param ts, time step for integration of \code{brownian.bridge.dyn()} object
#' @param ras, raster background size for \code{brownian.bridge.dyn()} object
#' @param le, location error value for \code{location.error=} used in a typical dbbmm object
#' @param lev, level of the UD contour as a vector, c(.50, .95) that the user is interested in.   Will work with #' multiple values (e.g., c(.50, .)) and will label those values
#' in the resultant shapefile.  Note that you have to put these values in as decimal (e.g., 0.50 for #' 50 percent).
#' @param crs, coordinate reference system for identifying where the polygons are located.  Uses standard \code{CRS} structure within quotes: CRS("+proj=longlat +datum=NAD83")
#' @param path, file path (e.g., "C:/") using standard R path nomenclature specifying the location that the #'output shapefiles are to be written.  Note that the call to \code{writeOGR} does not
#' accept the path to use the final '/' in the path nomenclature (e.g. "C:/TestLocation is fine as path is being pasted into the \code{dsn=} portion of \code{writeOGR} or if you assign path= to a
#' non-existant folder, it will create a folder of that name for you and store the output there..
#' @param name, file name, in quotation's, specifying the output files name which will be written to the #'\code{path} designated above.  This does not
#' require a file extension (e.g., 'shp').
#' @param ID, Unique identifier for each individual that is added to each created shapefile for linking to other #'data structures
#' @return Output is a set of shapefiles, written to the specified \code{path}
#' @export
#' @references Kranstauber, D., R. Kays, S. D. Lapoint, M. Wikelski, and K. Safi.  2012.  A dynamic Brownian #bridge movement model to estimate
#' utlization distributions for heterogenous animal movement.  Journal of Animal Ecology 81: 738-746. DOI: 10
#' .1111/j.1365-2656.2012.01955.x.
#' @references Byrne, M. E., J. C. McCoy, J. Hinton, M. J. Chamberlain, and B. A. Collier. 2014. Using dynamic #'brownian bridge movement modeling to measure temporal
#' patterns of habitat selection. Journal of Animal Ecology, 83: 1234-1243. DOI:  10.1111/1365-2656.12205
#' @author Bret A. Collier <bret@@lsu.edu>
#' @keywords moveud

move.forud=function(x, range.subset, ts, ras, le, lev, crs, path, name, ID)
{
  object<-x@DBMvar
  if(length(range.subset)<2)
    stop("\nrange.subset must have length >1\n")
  if(length(le)==1) location.error=rep(c(le), nrow(object))
  if(length(le)>1)  location.error=c(le)
  for(i in range.subset)
  {
    object@interest<-rep(F, nrow(object)); object@interest[i]<-T;times=object@timestamps[i];var=object@means[i];
    x.out <- brownian.bridge.dyn(object,raster=ras, time.step=ts, location.error=location.error)
    xx=raster2contour(x.out, level=lev)
    xx=spTransform(xx, CRS=CRS('+proj=longlat +datum=WGS84'))
    xx=SpatialLines2PolySet(xx)
    xx=PolySet2SpatialPolygons(xx)
    xx=as(xx, "SpatialPolygonsDataFrame")
    xx$levels=lev
    xx$times=times
    xx$stepvar=var
    xx$ID=ID
    writeOGR(xx, dsn=path, layer=paste(name, i, sep=""), driver="ESRI Shapefile")
  }
}

#----- Compute sedimentation from bathymetry data
# Parameters
#   secteur: sector to consider
compute_sedimentation <- function(sector, path = ".", path_sect = ".", ext = NULL, res = 5, fun = "min", pb_center = c("center","both","pb"), polygons = NULL, years = NULL, interpolate = "tin", ...)
# pb_center : 
{
    pb_center <- match.arg(pb_center)
    if (is.null(polygons)){
       allSects <- readOGR(dsn = sprintf("%s/Secteurs de sondage.shp",path_sect))
       allSects <- spTransform(allSects,CRS("+init=epsg:2950"))
       polygons <- allSects[grep(sector,allSects$name),]          
    }
    if (is.null(ext)) ext <- extent(polygons)
    #--- Find the sounding to keep (1st of spring and last of autumn)
    file_list <- list.files(path, pattern = tolower(sector))
    file_list <- file_list[grep("txt",file_list)]
    Dind <- duplicated(substr(file_list,1,12),fromLast=T)
    file_list <- file_list[!Dind]
    objnam <- substr(file_list,1,12)
    objnamDF <- data.frame(year = substr(objnam,1,2), sector = substr(objnam,3,5), 
      type = substr(objnam,6,6), rank = substr(objnam,7,7), day = substr(objnam,8,10), 
      stringsAsFactors=FALSE)
    if (is.null(years)) years <- unique(objnamDF$year)
    years <- formatC(years,width=2,flag="0")
    inYears <- objnamDF$year %in% years
    file_list <- file_list[inYears]
    objnam <- objnam[inYears]
    objnamDF <- objnamDF[inYears,]
    dates <- strptime(paste(objnamDF$year, objnamDF$day, sep="-"), "%y-%j")
    ##### We determine the first sounding of the season
    # For the first sounding of the season, we exclude type 4 (after dredging), type 8 (end of season), and type 9
    typesFirst <- !objnamDF$type %in% c(4, 8, 9)
    firstInds <- aggregate(as.numeric(objnamDF$day[typesFirst]), 
      by = list(year = objnamDF$year[typesFirst]), 
      FUN = min)
    # We group together first soundings since it can be spread on several days (based on numbers)
    firstInds <- lapply(seq_len(nrow(firstInds)), function(i) firstInds[i,])
    firstClusts <- lapply(firstInds,function(x){
        levnum <- objnam[which(objnamDF$year %in% x[1] & objnamDF$day %in% x[2])]
        levnum <- levnum[1]
        return(which(substr(objnam,1,7) == substr(levnum,1,7)))
    })
    # We determine the corresponding dat
    firstDates <- dates[unlist(lapply(firstClusts, function(x){x[order(objnamDF$day[x])][1]}))]
    ##### We apply the same process for the last sounding
    # We hereby exclude type 1 (after winter season) and 3 (before dredging)
    typesLast <- !objnamDF$type %in% c(1, 3, 9)
    lastInds <- aggregate(as.numeric(objnamDF$day[typesLast]), 
      by = list(year = objnamDF$year[typesLast]), FUN = max)
    lastInds <- lapply(seq_len(nrow(lastInds)), function(i) lastInds[i,])
    lastClusts <- lapply(lastInds,function(x){
        levnum <- objnam[which(objnamDF$year %in% x[1] & objnamDF$day %in% x[2])]
        if (length(levnum) > 1) levnum <- levnum[length(levnum)]
        return(which(substr(objnam,1,7) == substr(levnum,1,7)))
    })
    names(lastClusts) <- lastInds$year
    lastDates <- dates[unlist(lapply(lastClusts,
        function(x){x[order(objnamDF$day[x], decreasing=TRUE)][1]}))]
    ##### We read the identified soundings
    idFiles <- unique(c(unlist(firstClusts,use.names=F),unlist(lastClusts,use.names=F)))
    filesToRead <- file_list[idFiles]
    nfiles <- length(filesToRead)
    fun <- rep_len(fun,2)
    sectorRaster <- vector("list", nfiles)
    names(sectorRaster) <- substr(filesToRead,1,12)
    print(sprintf("sector %s: lecture de %s fichiers",sector,nfiles))
    pb <- txtProgressBar(min = 0, max = nfiles, style = 3)
    for (i in 1:nfiles){
        dat <- read.table(sprintf("%s/%s",path,filesToRead[i]))
        if (idFiles[i] %in% unlist(firstClusts, use.names=F)){
           funi <- fun[1]
        } else {
           funi <- fun[2]
        }
        # This custom function reads a bathymetric sounding in txt and transform it as a raster
        sectorRaster[[i]] <- bathyToRaster(dat, res = res, fun = funi, ext = ext)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    ##### Put first and last soundings together
    firstSounds <- stack(lapply(firstClusts,function(x){
        inds <- which(idFiles %in% x[order(objnamDF$day[x])])
        Reduce(cover,sectorRaster[inds]) # Paste together the rasters from the same sounding
    }))
    lastSounds <- stack(lapply(lastClusts,function(x){
        inds <- which(idFiles %in% x[order(objnamDF$day[x], decreasing=T)])   
        Reduce(cover,sectorRaster[inds]) # Paste together the rasters from the same sounding
    }))
    ##### Compute the raster difference
    # Years for which we have both an autumn and a spring sounding
    availYears <- intersect(firstDates$year, lastDates$year + 1)
    availYearsC <- formatC(availYears - 100, width = 2, flag = 0)
    # Initilization of rasterBrick
    sedimentRaster <- brick(nrows=nrow(sectorRaster[[1]]), ncols=ncol(sectorRaster[[1]]), 
      xmn=extent(sectorRaster[[1]])[1], xmx=extent(sectorRaster[[1]])[2], 
      ymn=extent(sectorRaster[[1]])[3], ymx=extent(sectorRaster[[1]])[4], 
      nl=length(years))
    datDF <- data.frame(Autumn = rep(as.POSIXct(NA), length(years)), 
      Spring = rep(as.POSIXct(NA),length(years)), 
      Diff = rep(NA,length(years)))
    if (length(availYears) > 0){
      diffRaster <- firstSounds[[which(firstDates$year %in% availYears)]] - 
        lastSounds[[which((lastDates$year + 1) %in% availYears)]] # The + 1 to match autumn with spring
      ndays <- round(julian(firstDates[firstDates$year %in% availYears]) - 
        julian(lastDates[(lastDates$year+1) %in% availYears]))
      for (i in seq_along(availYears)){
          sedimentRaster[[which(years == availYearsC[i])]] <- diffRaster[[i]]
          datDF[which(years == availYearsC[i]),]$Autumn <- lastDates[(lastDates$year+1) %in% availYears[i]]
          datDF[which(years == availYearsC[i]),]$Spring <- firstDates[firstDates$year %in% availYears[i]]
          datDF[which(years == availYearsC[i]),]$Diff <- ndays[i]
      }
      names(sedimentRaster) <- sprintf("A%s",years) 
      # We exclude sides
      poly_keep <- switch(pb_center,
          center = !grepl("[_(][SN][_)]",polygons$name),
          both = rep(TRUE,length(polygons)),
          pb = grepl("[_(][SN][_)]",polygons$name)
      )
      sedimentRaster <- if (sum(poly_keep) > 0) mask(sedimentRaster,polygons[poly_keep,]) else NA
    }  
    return(list(raster = sedimentRaster, dates = datDF))
}



#! Read bathymetric data
read_bathy <- function(path = ".", select = list(), funs = NULL, arguments = list(), meta.select = NULL)
# path : character string, the path indicating where the files are ;
# select : a list indicating the restrictions in the files read. For instance select  = list(sector = "d16") only reads the data from sector d16 ;
# funs : functions to apply to the basic bathymetry matrix. Useful to extract particular features ;
# arguments : optional arguments to be passed to the functions in funs. Each element of arguments must be named after a function in aggreg and is a list containing the values of arguments ;
# meta.select : the metadata fields to be returned. If NULL, no metadata is read.
{
    # Lowering case
    select <- lapply(select,tolower)
    names(select) <- sapply(names(select),tolower)
    pattern_components <- list(annee = "[[:digit:]][[:digit:]]", secteur = "[[:alpha:]][[:digit:]][[:digit:]]", type = "[[:digit:]]", rang = "[[:digit:]]", jour = "[[:digit:]][[:digit:]][[:digit:]]", navire = "[[:alnum:]]", num = "[[:digit:]]", ext = ".*.txt")
    # Checking that elements of select are legit
    check_names <- names(select) %in% names(pattern_components)
    if (any(!check_names)){
       warning(sprintf("Elements (%s) in select not used",paste(names(select)[!check_names],collapse=",")))
    }
    # Checking that elements are character
    is_char <- sapply(select,is.character)
    if (any(!is_char)){
       select[!is_char] <- lapply(select[!is_char],as.character)
       warning(sprintf("Elements (%s) coerced into character",paste(names(select)[!is_char],collapse=",")))
    }
    # creating regular expression 
    to_replace <- names(pattern_components) %in% names(select) 
    pattern_components[names(pattern_components)[to_replace]] <- select[names(pattern_components)[to_replace]]
    # checking if several choices are provided
    for (i in which(to_replace)){
        if (length(pattern_components[[i]]) > 1) pattern_components[[i]] <- sprintf("(%s)",paste(pattern_components[[i]],collapse="|"))
    }
    pattern <- paste(pattern_components,collapse="")    
    file_list <- list.files(path, pattern = pattern)
    nfiles <- length(file_list)
    print(sprintf("Reading %s soundings",nfiles)); flush.console()
    ids <- gsub(".txt","",file_list)
    nf <- length(funs)
    bathy <- vector("list",nfiles)
    names(bathy) <- ids
    nm <- length(meta.select)
    if (nm > 0){
       meta <- as.data.frame(matrix(NA, nfiles, nm))
       names(meta) <- meta.select
    } else {
       meta <- NULL
    }
    pb <- txtProgressBar(min = 0, max = nfiles, style = 3)
    for (i in 1:nfiles){
        dat <- read.table(sprintf("%s/%s",path,file_list[i]))
        if (nf > 0){
           funList <- vector("list",nf)
           names(funList) <- funs
           for (a in 1:nf){
               pars <- c(list(x = dat),arguments[[funs[a]]])
               funList[[a]] <- do.call(funs[a],pars)
           }
           if (nf == 1) funList <- funList[[1]]
           bathy[[ids[i]]] <- funList
        } else {
           bathy[[ids[i]]] <- dat
        } 
        #! Metadata
        if (nm > 0){
          meta_i <- xmlRoot(xmlParse(sprintf("%s/%s",path,sub("txt","object.xml",file_list[i]))))
          vals <- xmlApply(meta_i,xmlValue)
          names(vals) <- xmlApply(meta_i,xmlAttrs)
          foundAttr <- meta.select[meta.select %in% names(vals)]
          meta[i,foundAttr] <- unlist(vals[foundAttr])
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)    
    return(list(values = bathy, metadata = meta))
}

#! Count bathymetric data
count_bathy <- function(path = ".", by = c(1:7))
# path : character string, the path indicating where the files are ;
# by : same as the argument by in the function aggregate.
{
    file_list <- list.files(path, pattern = "[[:digit:]][[:digit:]][[:alpha:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:alnum:]][[:digit:]].txt")
    file_split <- id_components <- list(annee = c(1,2), secteur = c(3,5), type = c(6,6), rang = c(7,7), jour = c(8,10), navire = c(11,11), num = c(12,12))
    for (i in 1:length(file_split)) file_split[[i]] <- substr(file_list,id_components[[i]][1],id_components[[i]][2])
    file_split <- as.data.frame(file_split)
    n <- 1:length(file_list)
    count <- aggregate(n, by=file_split[,by,drop=F], length)
    return(count)
}

#! Transform bathymetric data into raster
bathyToRaster <- function(x, res, ext = extent(x), ...)
# x : reading from bathymetric data
# resolution : the resolution of the ratser (in meters)
# ... : arguments to pass to the rasterize function
{
    names(x) <- c("x","y","value")
    out <- raster(ext = ext)
    res(out) <- res
    out <- rasterize(x[,-3], out, x[,3], ...)
}

#! Locate secteur
locateSecteur <- function(object)
# n : number of locations to ping
# object : spatial object containing the limits of the sectors
{
    positions <- as.data.frame(locator(2))
    positions <- expand.grid(positions)[c(1:2,4:3),]
    rectangle <- SpatialPolygons(list(Polygons(list(Polygon(positions)),ID="rect")),proj4string=CRS("+init=epsg:2950"))
    covered <- intersect(rectangle,object)$name
    covered <- unique(substr(covered,1,3))
    return(covered)
}


#! Obtain sedimentation series



aggregate_sedimentation <- function(sedimentationRaster, increment, centre, fun = 'mean')
# sedimentationRaster : raster de sedimentations
# increment : en mètre la différence entre deux points de centre
# centre : la ligne de centre du chenal au long de laquelle les points sont tirés
{
  suppressWarnings(centre <- intersect(centre,as(extent(sedimentationRaster),"SpatialPolygons")))
  # Grille régulière de points sur la ligne de centre
  npts <- gLength(centre) / increment
  centrePts <- spsample(centre,n=npts,type="regular")
  npts <- length(centrePts)
  # Point le plus proche
  sedPts <- rasterToPoints(sedimentationRaster)
  near <- nncross(ppp(sedPts[,1],sedPts[,2],window=owin(extent(sedimentationRaster)[1:2],extent(sedimentationRaster)[3:4])),ppp(centrePts@coords[,1],centrePts@coords[,2],window=owin(extent(sedimentationRaster)[1:2],extent(sedimentationRaster)[3:4])),what="which")
  whichPt <- raster(ext = extent(sedimentationRaster))
  res(whichPt) <- res(sedimentationRaster)
  whichPt <- rasterize(sedPts[,1:2],whichPt,field=near)
  # agrégation spatiale
  rm(sedPts,near)
  final.data <- zonal(sedimentationRaster, whichPt, fun = fun)  
   #! Premiere version de la fonction
#  # Point le plus proche pour chaque pixel
#  distances <- vector("list",npts)
#  for (i in 1:npts) distances[[i]] <- distanceFromPoints(sedimentationRaster,centrePts[i])
#  whichPt <- calc(stack(distances),fun=which.min)
#  # Aggrégation pour chaque points
#  final.data <- matrix(NA,npts,dim(mergedRaster)[3])
#  colnames(final.data) <- names(sedimentationRaster)
#  for (i in 1:npts){
#      rasti <- mask(mergedRaster,whichPt==i,maskvalue=0)
#      final.data[i,] <- cellStats(rasti,stat="mean")
#  }
  return(list(series = final.data, pts = centrePts))
}

aggregate_sedimentation_perpendicular <- function(sedimentationRaster, increment, centre, secteurs, perp.increment = 20, fun = 'mean')
# sedimentationRaster : raster de sedimentations
# increment : en mètre la différence entre deux points de centre
# centre : la ligne de centre du chenal au long de laquelle les points sont tirés
{
  suppressWarnings(centre <- intersect(centre,as(extent(sedimentationRaster),"SpatialPolygons")))
  suppressWarnings(secteurs <- intersect(secteurs,as(extent(sedimentationRaster),"SpatialPolygons")))
  # Grille régulière de points sur la ligne de centre
  npts <- gLength(centre) / increment
  centrePts <- spsample(centre,n=npts,type="regular")
  npts <- length(centrePts)
  # Point le plus proche
  sedPts <- rasterToPoints(sedimentationRaster)
  spPts <- SpatialPoints(sedPts[,1:2], proj4string = centre@proj4string)
  distToCentre <- gDistance(spPts, centre, byid = TRUE) 
  distToSouth <- gDistance(spPts, secteurs[grep("\\(S\\)",secteurs$name),], byid = c(TRUE, TRUE))
  distToSouth <- apply(distToSouth, 2, min)
  distToNorth <- gDistance(spPts, secteurs[grep("\\(N\\)",secteurs$name),], byid = c(TRUE, TRUE))
  distToNorth <- apply(distToNorth, 2, min)
  distToCentre[distToNorth < distToSouth] <- -distToCentre[distToNorth < distToSouth]   
  near <- nncross(ppp(sedPts[,1],sedPts[,2],window=owin(extent(sedimentationRaster)[1:2],extent(sedimentationRaster)[3:4])),ppp(centrePts@coords[,1],centrePts@coords[,2],window=owin(extent(sedimentationRaster)[1:2],extent(sedimentationRaster)[3:4])),what="which")
  perp.dist <- as.numeric(as.character(cut(distToCentre, breaks = perp.increment, labels = 1:perp.increment)))
  near.perp <- 1000*near + perp.dist  
  whichPt <- raster(ext = extent(sedimentationRaster))
  res(whichPt) <- res(sedimentationRaster)
  whichPt <- rasterize(sedPts[,1:2],whichPt,field=near.perp)  
  # agrégation spatiale
  rm(sedPts)
  final.data <- zonal(sedimentationRaster, whichPt, fun = fun)
  nnear <- length(unique(near)) 
  curves <- lapply(1:nnear, function(x){ 
    temp <- final.data[floor(final.data[,1]/1000) == x,, drop = F]
    cur <- matrix(NA, perp.increment, ncol(final.data) - 1)
    inds <- temp[,1] %% 1000
    cur[inds,] <- temp[,-1]
    return(cur)
  })  
  return(list(curves = curves, pts = centrePts))
}


#! Test de corrélation par permutation
corperm <- function(x, y, N=1000, ...){  
    reps <- rep(NA,N)
    for (i in 1:N) reps[i] <- cor(sample(x),y,...)
    obs <- cor(x,y,...)
    p <- mean(abs(reps) > abs(obs))
    return(list(rho.obs = obs, p.value = p))
}

#! Plot des points
plot.pts <- function(x, values = NULL, borders = NULL, secteurs = NULL, increment = min(dist(x@coords)), asp.type = c("normal","horizontal","vertical"), pal = c("blue","white","red"), nc = 20, colorbar = TRUE, args.plot = list(), args.limits = list(pos = c(1,3), xpd = T, labels = c("Montréal","Trois-Rivières"), cex = 1.5), args.secteurs = list(pos = 1), args.colorbar = list(), args.borders = list())
{
    asp.type <- match.arg(asp.type)
    coords <- x@coords
    dc <- as.matrix(dist(coords))[1,]
    coords <- switch(asp.type,
        normal = coords,
        horizontal = data.frame(x = dc, y = rep(0, length(x))),
        vertical = data.frame(x = rep(0, length(x)), y = dc)
    )
    ncoords <- nrow(coords)
    if (!is.null(values)){
      colpal <- if (is.function(pal)) pal else colorRampPalette(pal,space="rgb")
      if (is.null(args.colorbar$zlim)) args.colorbar$zlim <- range(values, na.rm = T)
      colbr <- seq(args.colorbar$zlim[1], args.colorbar$zlim[2], 
        length.out = nc + 1)
      colvec <- colpal(nc)[cut(as.numeric(values), colbr)]
    } else {
      colvec <- pal[1]
    }
    args.plot$x <- coords
    args.plot$col <- colvec
    do.call(plot,args.plot)
    if (!is.null(borders)){
       args.borders$x <- borders
       args.borders$add <- TRUE
       do.call(plot,args.borders)
       do.call(points,args.plot[names(args.plot) %in% c("x","y",names(formals(plot.xy)))])
    }
    sc <- apply(coords,1,sum)
    inds <- c(which.min(sc),which.max(sc))
    args.limits$x <- coords[inds,]
    do.call(text,args.limits)
    if (!is.null(secteurs)){
       pt.sect <- over(x,secteurs)$name
       centr <- aggregate(coords,by = list(pt = substr(pt.sect,1,3)),mean)
       args.secteurs$x <- centr[,-1]
       args.secteurs$y <- centr[,1]
       do.call(text,args.secteurs)
    }
    if (colorbar){
       if (is.null(args.colorbar$zlim)) args.colorbar$zlim <- range(values)
       args.colorbar$col <- colpal(nc)
       do.call(image.scale,args.colorbar)
    }
}



#------- Draw a linear scale on the margin of the plot 
image.scale <- function(zlim, col = heat.colors(12), horiz = F, margin = T, zlab="", type.labels = c("continuous", "classes"), reverse = F, formatC.pars = list(), tck = 0.01, ticks.pars = list(), axis.pars = list(), label.pars = list())
# zlim: limits of the scale
# col: vector of colors of the scale. Also gives the number of breaks of the scale
# horiz: if TRUE the scale is drawn horizontally, if FALSE drawn vertically
# margin: if TRUE the scale is drawn in the margin of an existing plot, if false it is drawn in a new plot
# zlab: the scale label
# type.label: the type of label to draw
# reverse: if TRUE, the scale is reversed
# formatC.pars: parameters to format number in the scale
# tck: ticks length
# tick.pars: other parameters for ticks
# axis.pars: parameters for the scale legend text
# label.pars: parameters for the label
{
 type.labels <- match.arg(type.labels)
 if (!margin) plot.new()
 # side of the scale
 brpl_sides <- if (horiz) c(1,2) else c(2,1)
 brlim <- par("usr")[(brpl_sides[1] - 1) * 2 + 1:2]
 # Compute linear breaks and polygon limits
 breaks <- seq(brlim[1], brlim[2], length.out=(length(col)+1)) 
 poly <- vector(mode="list", length(col))
 for (i in seq(poly)){
     poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 # Prepare ticks
 tickl <- tck * min(par("usr")[4]-par("usr")[3],par("usr")[2]-par("usr")[1])
 formatC.pars$x <- switch(type.labels,
      continuous = c(zlim,mean(zlim)),
      classes = seq(zlim[1],zlim[2],length.out=length(poly))
 )
 labs <- do.call("formatC", formatC.pars)
 if (reverse){
    labs <- rev(labs)
    col <- rev(col)
 }
 # Plot the scale
 pllim <- par("usr")[(brpl_sides[2] - 1) * 2 + 1:2]
 if (margin){
    usrmar <- diff(pllim) * par("mai")[brpl_sides[2]+2] / par("pin")[brpl_sides[1]]
    pllim <- pllim[2] + usrmar * c(0.2,0.4)
 }
 poly.pars <- list(x = rep(pllim, each = 2), y = rep(pllim, each = 2), border = NA, xpd = NA) 
 for (i in seq(poly)){
     poly.pars[[brpl_sides[1]]] <- poly[[i]]
     poly.pars$col <- col[i]
     do.call("polygon", poly.pars) 
 }
 poly.pars[[brpl_sides[1]]] <- c(brlim[1],brlim[2],brlim[2],brlim[1])
 poly.pars$col <- NA
 poly.pars$border <- "black"
 do.call("polygon", poly.pars)
 # Add ticks
 axis.pars$labels <- labs
 axis.pars$pos <- brpl_sides[1] + 2
 axis.pars$xpd <- NA
 if (type.labels == "continuous"){
    ticks.pars$xpd <- NA
    ticks.pars[c("x0", "y0", "x1", "y1")] <- if (horiz) list(c(brlim,mean(brlim)),pllim[2],c(brlim,mean(brlim)),pllim[2]+tickl) else list(pllim[2],c(brlim,mean(brlim)),pllim[2]+tickl,c(brlim,mean(brlim)))
    do.call("segments", ticks.pars)
    axis.pars[c("x","y")] <- list(c(brlim,mean(brlim)),pllim[2]+tickl)[brpl_sides]    
 } else {
    axis.pars[c("x","y")] <- list(sapply(poly,function(x)mean(x[1:2])),pllim[2]+tickl)[brpl_sides]
 }
 do.call("text", axis.pars)
 # Add axis label
 lfun <- ifelse(horiz, "strheight", "strwidth")
 labwdth <- max(do.call(lfun, list(labs)))
 label.pars$xpd <- NA
 label.pars$label <- zlab
 if (horiz){
    label.pars[c("x","y","pos","offset")] <- list(mean(brlim),pllim[2]+tickl+labwdth,3,1.5)
 } else {
    labpos <- grconvertX(pllim[2]+tickl+labwdth,from="user",to="inches")
    laboff <- strheight(zlab,units="inches")
    label.pars[c("x","y","srt")] <- list(grconvertX(labpos+2*laboff,from="inches",to="user"),mean(brlim),270)
 }
 do.call("text", label.pars)
}
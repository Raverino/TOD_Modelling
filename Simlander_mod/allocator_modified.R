#begin script
#set workspace and libraries
#setwd("C:/simlandertfm/SIMLANDERICHARD")
#require(raster)
#require(rgdal)
#require(gWidgets2tcltk)

#library(sf)
#library(igraph)

library(terra)
library(landscapemetrics)

####ALL LIBRARIES PREVIOUSLY LOADED IN SIMCONTROL. IF NOT (I.E. IN CASE OF NOT BEING ABLE TO USE GUI), UNCOMMENT TO LOAD LIBRARIES SEPARATELY IN allocator.R
#set.seed(102018)   ##SET SEED NOW MANAGED FROM SIMCONTROL DIALOG. 

##get startyear and endyear (previously set from simcontrol dialog)

####FUNCTION MAPCOMPARE
mapcompare <- function(Tn,r1,runname) { ##function to compare Tn with reference map for interation i
	
	statline <- list()  ##create the list and datframe to receive the stats (only one line)
	statline <- as.data.frame(statline)  ##create the list and datframe to receive the stats (only one line)
	
	s1 <- Tn
	s1test <- s1
	#s1test[is.na(s1test[])] <- 0
	s1[s1%in%c(1)] <- 999   #urban land in the simulation is 999
	s1r1 <- s1+r1
	xvals <- as.data.frame(values(s1r1))
	colnames(xvals) <- "values"
	MISSES <- sum(xvals$values == 1,na.rm=T)
	FALSE_ALARMS <- sum(xvals$values == 999,na.rm=T)
	HITS <- sum(xvals$values == 1000,na.rm=T)   #value 1 in ref map, value 999 in simulated map
	recall <- HITS/(HITS+MISSES)
	
	####CALCULATE CLASS LEVEL STATS simulations##########landscapemetrics package
	ED <- as.numeric(lsm_l_ed(s1)[6]) #Edge Density
	#EDna <- as.numeric(lsm_l_ed(s1test)[6]) #Edge Density
	FD <- as.numeric(lsm_c_frac_mn(s1)[2,6]) #Mean Fractal Dimension Index
	#FDna <- as.numeric(lsm_c_frac_mn(s1test)[2,6]) #Mean Fractal Dimension Index
	TCA <- as.numeric(lsm_c_tca(s1)[2,6]) #Total Core Area
	#TCAna <- as.numeric(lsm_c_tca(s1test)[2,6]) #Total Core Area
	CL <- as.numeric(lsm_c_clumpy(s1)[2,6]) #Clumpiness Index
	#CLna <- as.numeric(lsm_c_clumpy(s1test)[2,6]) #Clumpiness Index
	PD <- as.numeric(lsm_c_pd(s1)[2,6]) #Patch Density
	#PDna <- as.numeric(lsm_c_pd(s1test)[2,6]) #Patch Density
	####end CALCULATE CLASS LEVEL STATS simulations##########landscapemetrics package
	
	##Recall = Hits / (Hits + Misses) 
	statline[1,1] <- runname	##name of the map from the list
	statline[1,2] <- HITS
	statline[1,3] <- MISSES
	statline[1,4] <- FALSE_ALARMS
	statline[1,5] <- recall
	statline[1,6] <- ED
	#statline[1,7] <- EDna
	statline[1,7] <- FD
	#statline[1,9] <- FDna
	statline[1,8] <- TCA
	#statline[1,11] <- TCAna
	statline[1,9] <- CL
	#statline[1,13] <- CLna
	statline[1,10] <- PD
	#statline[1,15] <- PDna
	#colnames(statline) <- c("mapname","hits","misses","false alarms","recall (H/(H+M)")
	return(statline)
	#list(statline=statline,refline=refline)
}


####END FUNCTION MAPCOMPARE

if (file.exists("INPUTS/startyear.txt")) {
startyear <- as.numeric(readLines("INPUTS/startyear.txt"))
} else {
tkmessageBox(message = "please run simcontrol to set startyear")
}

if (file.exists("INPUTS/endyear.txt")) {
endyear <- as.numeric(readLines("INPUTS/endyear.txt"))
} else {
tkmessageBox(message = "please run simcontrol to set startyear")
}

if (file.exists("INPUTS/modeltype.txt")) {
mtype <- readLines("INPUTS/modeltype.txt")
} else {
tkmessageBox(message = "please run simcontrol to set model type")
}

if (file.exists("INPUTS/rantype.txt")) {
rantype <- readLines("INPUTS/rantype.txt")
} else {
tkmessageBox(message = "please run simcontrol to set random seed type")
}

if (file.exists("INPUTS/auto.txt")) {
ardswitch <- as.numeric(readLines("INPUTS/auto.txt"))
} else {
tkmessageBox(message = "please run simcontrol to define autocalibration settings")
}

if (file.exists("INPUTS/numruns.txt")) {
numruns <- as.numeric(readLines("INPUTS/numruns.txt"))
} else {
tkmessageBox(message = "please run simcontrol to set number of runs")
}

if (file.exists("loopval.txt")) {
loopval <- as.numeric(readLines("loopval.txt"))
} else {
tkmessageBox(message = "please run simcontrol to create loopval.txt")
}

#get startmap
if (file.exists("INPUTS/inmap.txt")) {
	lu1path  <- readLines("INPUTS/inmap.txt")
	print(paste0("got initial map for year ",startyear," from file ",lu1path))
	lu1 <- rast(lu1path)
	crs(lu1) <- "epsg:25830"
	v <- values(lu1)
	vf <- as.data.frame(v)
	vf <- na.omit(vf)
	##urbdemand <- sum(vf$v)   #for the raster package
	urbdemand <- sum(vf)
} else {
tkmessageBox(message = "please run simcontrol to define start map")
}	

if (file.exists("INPUTS/outmap.txt")) {
	lu2path  <- readLines("INPUTS/outmap.txt")
	print(paste0("got final map for year ",endyear," from file ",lu2path))
	lu2 <- rast(lu2path)
	crs(lu2) <- "epsg:25830"
	r1 <- lu2   ###if second reference map has been included, statsvals will compare against it
	#lu1[lu1%in%c(1)] <- 0
	#lu1[lu1%in%c(11)] <- NA
	#lu1[lu1 > 0] <- 1
	v2 <- values(lu2)
	vf2 <- as.data.frame(v2)
	vf2 <- na.omit(vf2)
	#lu2cells <- sum(vf2$v2)  #for the raster package
	lu2cells <- sum(vf2)
} else {
r1 <- lu1 ###if no second reference map has been included, statsvals will compare against lu1
}


if (file.exists("INPUTS/statsvals.txt")) {
	library(landscapemetrics)
	statsvals <- as.numeric(readLines("INPUTS/statsvals.txt"))
	##3source("INPUTS/mapcompare.R")  #load function to compare maps on the fly
	statlist <- list()  ##create the lists and datframe to receive the stats
	statlist <- as.data.frame(statlist)  ##create the lists and datframe to receive the stats
	refline <- list()  ##create the list and datframe to receive the stats (only one line)
	refline <- as.data.frame(refline)  ##create the list and datframe to receive the stats (only one line)
	####CALCULATE CLASS LEVEL STATS reference map##########landscapemetrics package
	EDr <- as.numeric(lsm_l_ed(r1)[6]) #Edge Density
	FDr <- as.numeric(lsm_c_frac_mn(r1)[2,6]) #Mean Fractal Dimension Index
	TCAr <- as.numeric(lsm_c_tca(r1)[2,6]) #Total Core Area
	CLr <- as.numeric(lsm_c_clumpy(r1)[2,6]) #Clumpiness Index
	PDr <- as.numeric(lsm_c_pd(r1)[2,6]) #Patch Density
	####end CALCULATE CLASS LEVEL STATS reference map##########landscapemetrics package
	####reference map comparisons
	refline[1,1] <- paste0("Reference Map r1 = ",names(r1))	##name of the map from the list
	refline[1,2] <- NA
	refline[1,3] <- NA
	refline[1,4] <- NA
	refline[1,5] <- NA
	refline[1,6] <- EDr
	refline[1,7] <- FDr
	refline[1,8] <- TCAr
	refline[1,9] <- CLr
	refline[1,10] <- PDr
	######end reference map comparisons 
	
} else {
	tkmessageBox(message = "map statistics not set, if you need them, run again & choose 'output map comparison statistics' in the simcontrol dialog")
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#import and calculate accessibility
if (file.exists("INPUTS/accmap.txt")) {
	accpath  <- readLines("INPUTS/accmap.txt")
	distroad <- rast(accpath)  ####add in the roadnetwork 
	crs(distroad) <- "epsg:25830"
	#distroad <- resample(distroad, lu1, method ="bilinear")  ##its a different size, hence resample
	#writeRaster(distroad,"RoadDistances.tif",format="GTiff",overwrite=T)
	distroad_null <- mask(distroad,lu1)
	#begin calculate accessibility with White's accessibility equation
	accurb_road <- (1+(distroad/1))^-1
	
	# Add the new transit network accessibility map
	if (file.exists("INPUTS/transit_accmap.txt")) {
	  transit_accpath <- readLines("INPUTS/transit_accmap.txt")
	  disttransit <- rast(transit_accpath)  # Transit network
	  crs(disttransit) <- "epsg:25830"
	  disttransit_null <- mask(disttransit, lu1)
	  # Calculate transit accessibility with White's accessibility equation
	  accurb_transit <- (1+(disttransit/1))^-1
	
	
	  # Combine the two accessibility maps with 50% weight each
	  model_accessibility <- (0.2* accurb_road) + (0.8* accurb_transit)
	} else {
	  # If transit map doesn't exist, use only road accessibility
	  print("Transit network map not found, using only road accessibility")
	  model_accessibility <- accurb_road
	}
} else {
  tkmessageBox(message = "please run simcontrol to set distance to infrastructures map")
}


	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#import and calculate suitability
if (file.exists("INPUTS/suitmap.txt")) {
	suitpath  <- readLines("INPUTS/suitmap.txt")
	slope <- rast(suitpath)  ####add in the slope map 
	crs(slope) <- "epsg:25830"
	m <- c(-Inf,1,1,
	1,2,0.9,
	2,3,0.8,
	3,4,0.5,
	4,Inf,0.1)
	rclmat <- matrix(m, ncol=3, byrow=TRUE)
	recslope <- classify(slope, rclmat, include.lowest=TRUE)
	##recslope <- reclassify(slope,c(-Inf,1,1,
	##1,2,0.9,
	##2,3,0.8,
	##3,4,0.5,
	##4,Inf,0.1))
	
	recslope <- resample(recslope, lu1, method ="bilinear")
	slope_null <- mask(recslope,lu1)
	model_suitability <- slope_null						
	#end calculate slope
} else {
	tkmessageBox(message = "please run simcontrol to create a slope or other suitability map")
}


# Add zoning constraint after loading suitability 
# Zoning map should be:
# 1 = development allowed; 0 = development prohibited; NA = no zoning restrictions
# Don't add a zoningmap.txt if you don't want to add zoning, it will be identified as "non-existent" and run by itself
if (file.exists("INPUTS/zoningmap.txt")) {
  zoningpath <- readLines("INPUTS/zoningmap.txt")
  zoning <- rast(zoningpath)
  crs(zoning) <- "epsg:25830"
  zoning <- resample(zoning, lu1, method="near")  # Use "near" for categorical data
  
  model_zoning <- zoning
  model_zoning[is.na(model_zoning)] <- 1  # Treat NA as allowed
  
  # Making sure all urban cells are 0 or 1
  existing_urban <- lu1
  existing_urban[existing_urban != 1] <- 0  # Convert non-urban to 0
  existing_urban[existing_urban == 1] <- 1  # Keep urban as 1
  
  # Combine: Allow development where either zoning permits OR urban already exists
  model_zoning <- max(model_zoning, existing_urban, na.rm=TRUE)
  
  print("Zoning constraint created - existing urban areas preserved")
  zoning_exists <- TRUE
  
} else {
  print("No zoning map found - no zoning constraints applied")
  zoning_exists <- FALSE
}

	#get final land use demands
if (file.exists("INPUTS/demand.txt")) {
	finaldemand <- as.numeric(readLines("INPUTS/demand.txt"))
	} else {
		if (exists("lu2cells")) {
			tkmessageBox(message = "no specified final demand, using lu2")
			finaldemand <- lu2cells
		} else {
			tkmessageBox(message = "no specified final demand, please run simcontrol to create demand.txt or import map for endyear")
		}
	} 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##set seed according to whether we want the same or new random seed each time
eval(parse(text=rantype))  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
Tn <- lu1   ###important to define Tn before every simulation

if (ardswitch ==1) {
	####ARD matrix generator 
	####start the ARD matrix generator

	################################################################################
	#################################### Matrix ####################################
	################################################################################

	#read in the value of nw as created by the script file, so R doesn't run out of memory and crash
	nw <- 10   #set the number of random matrices to be generated
	
	size <- 11 # what is the size you want
	#building up random Pythagorian matrices

	randomPythagorianMatrix <- function(n, x, interpolation="smooth") {
  	seed  <- runif(n, 5, 50)
  	drops <- mapply(runif, 1, 1, seed)
  	Map(getPythagorianMatrix, x, seed, drops, interpolation)
	}

	getPythagorianMatrix <- function(x, mid, drop, interpolation="smooth") {
	if(x %% 2 == 0 | x < 0) stop("x must be an odd positive number")
  
	choices <- c("smooth", "linear")
	interpolation <- choices[pmatch(interpolation, choices, duplicates.ok=FALSE)]
  
	dists <- outer(abs(1:x - ceiling(x/2)), abs(1:x - ceiling(x/2)), function(x,y) sqrt(x^2+y^2))
	if(interpolation=="smooth") {
  	mat <- (1/drop) ^ dists * mid
  	} else {
  	mat <- matrix(approx(x=0:x, y=0.1^(0:x)*50, xout=dists)$y, ncol=x, nrow=x)
  	}
	return(mat)
	}


	############################## GENERATE MATRICES ############################################################# RANDOM MATRICES ################################
	mats <- randomPythagorianMatrix(nw, size, interpolation="smooth")
	
	#@@@cycle through a series a nhood rules selected from the matrices generated at the start of this script
	for (m in 1:length(mats)) {
	w.mrx <- mats[[m]]#mats[[nw]][c(4:8),c(4:8)]  
	
		for (j in 2:6)	{ 
			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			#set weights with nhood matrix
			#w2 <- mats[[m]]#mats[[nw]][c(4:8),c(4:8)] 					
			print (paste("simulation starting at ", Sys.time(), sep=""))
			total <- (endyear)-(startyear)  #set the total number of timesteps from endyear and startyear (set in simcontrol.R)

			##set annual demand (andem) using urbdemand - number of urban cells in startmap
			##and finademand, total final land claims required set in simcontrol.R
			andem <- (finaldemand-urbdemand)/total

			###start the timer....
			print("set demand, starting timer..")
			ptm <- proc.time()
			#end preparations
			#-----------------------
			pb <- txtProgressBar(min = 0, max = total, style = 3)
			
			for (i in 1:total){    #begin running simulation
				print (paste("simulation starting at ", Sys.time(), sep=""))
				#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@CURRENT DEMAND	
				##first set the amount of cells to be allocated in timestep i
				urbdemand1 <- urbdemand + andem*i 
	
				#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@NEIGHBOURHOOD BLOCK	
				#the calculate the neighbourhood block using the latest map Tn and the neighbiourhood weights  set above
				#begin neighbourhood block
      			w.s.rc <- (j-1) #first row/column of matrix to be considered for subsetting
      			w.e.rc <- 11-(j-2) #last row/column of matrix to be considered for subsetting
      
      			w <- w.mrx[c(w.s.rc:w.e.rc),c(w.s.rc:w.e.rc)]   
      			write.table(w, paste0("OUTPUTS/nhood_",m,"_r",7-j,".txt"),append = FALSE, col.names = NA,sep=",") 
      				
      			##plot the rules##########
      			matdims <- (length(w[,1])-1)/2
      			x <- 1:matdims * xres(lu1) ##the cell resolution
      			x0 <- 0
      			xax <- as.data.frame(c(x0,x))
      			x <- as.numeric(xax[,1])
      			centre <- w[matdims+1,matdims+1]
					for (a in 0:(matdims-1)) {
      					m1 <- w[matdims+1,matdims-a]
      					assign(paste0("ndist",a),m1)
      					}
      			dlist <- ls(pattern="ndist")
      			rules <- as.data.frame(c(centre,mget(dlist)))
      			rules <- as.numeric(rules[1,])
      			glmax <- max(unlist(lapply(mats,FUN=max)))  #get the global max to set yaxis height
      			png(filename = (paste0("OUTPUTS/nrules_",m,"_r",7-j,".png")),width = 600, height = 600, bg="white") #plot the rules
      			#plot(x,rules, type="o", col="blue",main = (paste0("nrules_",m,"_r",7-j)), xlab= "Distance (metres)", ylab= "influence score")
      			plot(x,rules, type="o", col="blue",main = (paste0("OUTPUTS/nrules_",m,"_r",7-j)), xlab= "Distance (metres)", ylab= "influence score",ylim=c(0,glmax))
      			dev.off()
      			rm(list=ls(pattern="ndist"))
      			##end plot the rules##########
      			#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
				#begin neighbourhood block 
				#n <- raster::focal(Tn, w=w2, na.rm=TRUE)   
				#nhood <- raster::cover(n, Tn)
				#model_nhood <- nhood
				
				n <- focal(Tn, w=w2, na.rm=TRUE)   
				nhood <- cover(n, Tn)
				model_nhood <- nhood
				#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
				#end neighbourhood block    
		
				#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@RANDOM BLOCK
				#begin random block    
				x <-runif(ncell(Tn)) #corrected to calculate directly from number of cells in the map
				weibull <- 1+(-log(1-(x)))*exp(1/2) ##set random factor here, e.g. exp(0.9)
				funselect <- function(x) { x[x!=0] <- NA; return(x) } #extract only vacant areas (value 0) so that existing functional land uses are not randomized. Set everything else to NA
				vacants <- calc(Tn, funselect)  
				#lu1 <- rast(lu1)
				#vacants <- app(lu1, funselect) 
				random <- lu1 #just copying lu1 as a skeleton map
				values(random)<-weibull 
				model_random <- mask(random, vacants) #generating a mask to apply NAs to the random layer. 
				model_random <- cover(model_random,Tn) #cover to fill the NAs back in with the original values from T1.
				#end random block   
	
				#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@TRANSITION POTENTIAL CALCULATION, ACCORDING TO MODEL TYPE
				CS <- lu1   #the CS allows us to give a strong weight to existing land use at time "startyear"
	
				##normalise model blocks, used by all models except TM_P
				maxn <- global(model_nhood, fun="max",na.rm=TRUE)
				minn <- global(model_nhood, fun="min",na.rm=TRUE)
				TPn <- (model_nhood - minn[1,]) / (maxn[1,] - minn[1,])

				maxr <- global(model_random, fun="max",na.rm=TRUE)
				minr <- global(model_random, fun="min",na.rm=TRUE)
				TPr <- (model_random - minr[1,]) / (maxr[1,] - minr[1,])
	
				maxa <- global(model_accessibility, fun="max",na.rm=TRUE)
				mina <- global(model_accessibility, fun="min",na.rm=TRUE)
				TPa <- (model_accessibility - mina[1,]) / (maxa[1,] - mina[1,])
				
				maxs <- global(model_suitability, fun="max",na.rm=TRUE)
				mins <- global(model_suitability, fun="min",na.rm=TRUE)
				TPs <- (model_suitability - mins[1,]) / (maxs[1,] - mins[1,])
				##end normalise model blocks,
				
				if (mtype == "PM") {
					print ("running the parsomonious model (inertia, nhood, random")
					model_TTP <- CS + TPn + TPr  #The Parsiomonius Model
					#end transition potential calculation  			
					########################################
					}
		
				else if (mtype == "AM") {
					print ("running the accessiblity model (inertia, nhood, random, accessibility")
					model_TTP <- CS + TPn + TPr + TPa #The Accessibility model
					#end transition potential calculation  			
					########################################
					}
	
				else if (mtype == "SM") {
					print ("running the suitability model (inertia, nhood, random, suitability")
					model_TTP <- CS + TPn + TPr + TPs #The Suitability Model
					#end transition potential calculation  			
					########################################
					}
	
				else if (mtype == "TM_S") {
					print ("running the typical model (inertia, nhood, random, accessibility suitability -  summed" )
					model_TTP <- CS + TPn + TPr + TPa + TPs  #The Typical Model 
					#end transition potential calculation  			
					########################################
					}

				else if (mtype == "TM_P") {
					print ("running the typical model (inertia, nhood, random, accessibility suitability -  multiplied" )
					model_TTP <- (model_accessibility+1)*(model_suitability+1)*(model_nhood+1)*(model_random)

				} else {
				print("model type not found, run simcontrol to set")
				}
	
				TTP <- mask(model_TTP,lu1)  #make sure we exclude the NAs 
				
				#TTP <- raster(TTP)
				#Tn <- raster(Tn)
				#lu1 <- raster(lu1)
	
    				########################################
   				# Begin Land use allocation #
    				########################################    ########################################		
			
   				x <- as.matrix(TTP,wide=T)  ##wide=T returns a matrix like 		that produced by the raster package
				#x <- as.matrix(TTP)  ##wide=T returns a matrix like 		that produced by the raster package
    			n <- urbdemand1  #the demand calculated as above
    			x2 <-  sort(-x, partial = n)
    			x2h <-  -sort(x2[1:n])
    			ix <- which(x %in% x2h)
    			rowsTn <- nrow(Tn)
    			colsTn <- ncol(Tn)
    			ro <- ix %% rowsTn
    			ro <- ifelse(ro == 0L, rowsTn, ro)
    			co <- 1 + ix %/% rowsTn
    			x3 <- x[ix]
    			d <- data.frame(row = ro, col = co, x = x3)
    			result <- d[rev(order(d$x)), ]
    				#-----------------------
				#test for duplicates that inflate the number of cells allocated
				difftrans <- (length(result$x)-n)
					if (difftrans > 0) {
					result2 <- head(result,-difftrans) #remove the duplicates from the end of the file (the weakest candidate cells)
					result <- result2
				}
				#turn selected n values in the dataframe into a matrix and then to a raster
				x.mat <- matrix(0, rowsTn, colsTn) #create a matrix with the right number of rows and columns and fill with 0 values
				#x.mat[cbind(result$row,result$col)] <-result$x   #works just fine. 
				x.mat[cbind(result$row,result$col)] <-1   #but actually we want values to be 1 (urban land), not the TP value. 
				#which(!is.na(x.mat), arr.ind =TRUE)  #pick out the non-NA values. Useful for testing that this has actually worked. 
				#r <- rast(x.mat)
				r <- rast(x.mat)
				#extent(r) <- extent(Tn)
				ext(r) <- ext(Tn)
				crs(r) <- crs(Tn)
				newdata <- r
				newdata <- mask(newdata, lu1) #generating a mask to apply NAs to the final map layer. 
				Tn <- newdata  #directly allocated all the cells at their most favourable locations according to TTP
				runname <- paste0(mtype,"_lu_", (startyear + i),"_run_",loopval,"_of_",numruns,"_m",m,"_r",7-j)
				filepng <- paste0("OUTPUTS/",runname,".png")
				#filenasc <- paste("lu", (startyear + i), ".asc", sep="")
				#LUout <- writeRaster(T1, filename=(filenasc), format="ascii", overwrite=TRUE)
				#filen <- paste(filen,".png", sep="")
				###@@@@@@@@@######png(filename = filepng,width = 1200, height = 1200, bg="white") #Plot each land use map to be able to make an animation.
				#plot(T1)
				#plot(T1,breaks=breakpoints,col=colors,main=(startyear+ii))
				####@@@@@@@#######plot(Tn,main=(startyear+i))
				####@@@@@@@#######dev.off()
		
				#removeTmpFiles(h=5)
				Sys.sleep(0.1)
				setTxtProgressBar(pb, i)
				print (paste("FINISHED: landuse simulation for ", (startyear + i), sep=""))
				writeRaster(Tn,paste0("OUTPUTS/",runname,"_",startyear+i,".tif"),overwrite=T)
				print(proc.time() - ptm)
			}
			close(pb)
			writeRaster(Tn,paste0("OUTPUTS/",runname,".tif"),overwrite=T)
				
			#######BLOCK TO CALCULATE MAP STATS ON THE FLY
				
			if (statsvals > 0) {
				print("creating map statistics table and running map statistics calculation...")
				statline <- mapcompare(Tn,r1,runname)
				statlist <- rbind(statlist,statline)   #add to the statslist file, initially empty, then updated as each map is produced
				##if (length(statlist) < 10)  {
				##	print("ERROR! can't bind statline to statlist")
				##	}
				} else {
				print("map statistics not required...")
			}
			print("map statistics calculation complete.")
			#######END BLOCK TO CALCULATE MAP STATS ON THE FLY

			sim <- Tn  
			
			#sss <- c(lu1,sim)
			#name1 <- paste("initial land use ",startyear,sep="")
			#name1 <- paste("land use in ",startyear,sep="")
			#name2 <- paste("simulated land use ",endyear,sep="")
			#names(sss) <- c(name1,name2) 
			#print(plot(sss)) #without the print command, the output will not appear on the screen
			print (paste("simulation complete at ", Sys.time(), sep=""))
			#print(plot(sss,breaks=breakpoints,col=colors))   #same as above, but just in case you have colors and breakpoints defined
			#save.image(paste0("sim_",loopval,".RData"))  #in case you want the workspace objects from each simulation
	
    		#naming variables
    		assign(paste("simT1_",j, sep=""), sim)
    		assign(paste("simT2_",j, sep=""), sim)
		}#end loop for rules sets j
	}#end loop for matrices m
	png(filename = paste("OUTPUTS/lu_",startyear,".png",sep=""),width = 1200, height = 1200, bg="white") #Plot initial map also to animations file
    plot(lu1, main = paste("initial land use map - ",startyear,sep=""))  #plot the start map  ## RH MODIFICATION (correcting old error by RH)
    dev.off()
	
} else {	
	
m <- 1 ##for the file output names
j <- 1
	
w2 <- matrix(c(0.000,0.000,0.500,0.000,0.000,
0.000,3.136,5.000,3.136,0.000,
0.500,5.000,50.000,5.000,0.500,
0.000,3.136,5.000,3.136,0.000,
0.000,0.000,0.500,0.000,0.000), 
            nr=5,nc=5)
            
#w2 <- w2*1000

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@			
				
print (paste("simulation starting at ", Sys.time(), sep=""))
total <- (endyear)-(startyear)  #set the total number of timesteps from endyear and startyear (set in simcontrol.R)

##set annual demand (andem) using urbdemand - number of urban cells in startmap
##and finademand, total final land claims required set in simcontrol.R
andem <- (finaldemand-urbdemand)/total

###start the timer....
print("set demand, starting timer..")
ptm <- proc.time()
#end preparations
#-----------------------
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (i in 1:total){    #begin running simulation
	print (paste("simulation starting at ", Sys.time(), sep=""))
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@CURRENT DEMAND	
	##first set the amount of cells to be allocated in timestep i
	urbdemand1 <- urbdemand + andem*i 
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@NEIGHBOURHOOD BLOCK	
	#the calculate the neighbourhood block using the latest map Tn and the neighbiourhood weights  set above   
	#n <- raster::focal(Tn, w=w2, na.rm=TRUE)   
	#nhood <- raster::cover(n, Tn)
	#model_nhood <- nhood
	
	n <- focal(Tn, w=w2, na.rm=TRUE)   
	nhood <- cover(n, Tn)
	model_nhood <- nhood
	
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
	#end neighbourhood block    
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@RANDOM BLOCK
	#begin random block    
	
	#set.seed(16785)
	x <-runif(ncell(Tn))
	weibull <- 1+((-log(1-(x)))*exp(1/2))
	funselect <- function(x) { x[x!=0] <- NA; return(x) } 
	vacants <- app(Tn, funselect) 
	random <- lu1
	values(random) <- weibull 
	model_random <- mask(random, vacants)
	model_random <- cover(model_random,Tn) 
	model_random <- model_random
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@TRANSITION POTENTIAL CALCULATION, ACCORDING TO MODEL TYPE
	CS <- lu1   #the CS allows us to give a strong weight to existing land use at time "startyear"
	
	##normalise model blocks, used by all models except TM_P
	maxn <- global(model_nhood, fun="max",na.rm=TRUE)
	minn <- global(model_nhood, fun="min",na.rm=TRUE)
	TPn <- (model_nhood - minn[1,]) / (maxn[1,] - minn[1,])

	maxr <- global(model_random, fun="max",na.rm=TRUE)
	minr <- global(model_random, fun="min",na.rm=TRUE)
	TPr <- (model_random - minr[1,]) / (maxr[1,] - minr[1,])
	
	maxa <- global(model_accessibility, fun="max",na.rm=TRUE)
	mina <- global(model_accessibility, fun="min",na.rm=TRUE)
	TPa <- (model_accessibility - mina[1,]) / (maxa[1,] - mina[1,])
				
	maxs <- global(model_suitability, fun="max",na.rm=TRUE)
	mins <- global(model_suitability, fun="min",na.rm=TRUE)
	TPs <- (model_suitability - mins[1,]) / (maxs[1,] - mins[1,])
	##end normalise model blocks,
	
	if (mtype == "PM") {
		model_TTP <- CS + TPn + TPr  #The Parsiomonius Model
		#end transition potential calculation  			
		########################################
		}
		
	else if (mtype == "AM") {
		model_TTP <- CS + TPn + TPr + TPa #The Accessibility model
	#end transition potential calculation  			
		########################################
		}
	
	else if (mtype == "SM") {
		model_TTP <- CS + TPn + TPr + TPs #The Suitability Model
	#end transition potential calculation  			
		########################################
		}
	
	else if (mtype == "TM_S") {
		model_TTP <- CS + TPn + TPr + TPa + TPs  #The Typical Model 
	#end transition potential calculation  			
		########################################
		}

	else if (mtype == "TM_P") {
		model_TTP <- (model_accessibility+1)*(model_suitability+1)*(model_nhood+1)*(model_random)
	#the typical SIMLANDER model   			
	#end transition potential calculation  			
		########################################

	} else {
	print("model type not found, run simcontrol to set")
	}
	
	TTP <- mask(model_TTP,lu1)  #make sure we exclude the NAs
	
	if (zoning_exists) {
	  TTP <- TTP * model_zoning    # Apply zoning constraint
	  print("Zoning constraints applied to transition potential")
	} else {
	  print("Running simulation without zoning constraints")
	}
	##TTP <- raster(TTP)
	##Tn <- raster(Tn)
	##lu1 <- raster(lu1)
	
	########################################
   	# Begin Land use allocation terra #
    ########################################   
	########################################		
	#TTP <- rast(TTP)
	#Tn <- rast(Tn)
   	x <- as.matrix(TTP,wide=T)  ##wide=T returns a matrix like 		that produced by the raster package
	#x <- as.matrix(TTP)  ##wide=T returns a matrix like 		that produced by the raster package
    n <- urbdemand1  #the demand calculated as above
    x2 <-  sort(-x, partial = n)
    x2h <-  -sort(x2[1:n])
    ix <- which(x %in% x2h)
    rowsTn <- nrow(Tn)
    colsTn <- ncol(Tn)
    ro <- ix %% rowsTn
    ro <- ifelse(ro == 0L, rowsTn, ro)
    co <- 1 + ix %/% rowsTn
    x3 <- x[ix]
    d <- data.frame(row = ro, col = co, x = x3)
    result <- d[rev(order(d$x)), ]
    #-----------------------
	#test for duplicates that inflate the number of cells allocated
	difftrans <- (length(result$x)-n)
	if (difftrans > 0) {
		result2 <- head(result,-difftrans) #remove the duplicates from the end of the file (the weakest candidate cells)
		result <- result2
		}
	#turn selected n values in the dataframe into a matrix and then to a raster
	x.mat <- matrix(0, rowsTn, colsTn) #create a matrix with the right number of rows and columns and fill with 0 values
				#x.mat[cbind(result$row,result$col)] <-result$x   #works just fine. 
	x.mat[cbind(result$row,result$col)] <-1   #but actually we want values to be 1 (urban land), not the TP value. 
	#which(!is.na(x.mat), arr.ind =TRUE)  #pick out the non-NA values. Useful for testing that this has actually worked. 
	r <- rast(x.mat)
	#r <- raster(x.mat)
	#extent(r) <- extent(Tn)
	ext(r) <- ext(Tn)
	crs(r) <- crs(Tn)
	newdata <- r
	newdata <- terra::mask(newdata, lu1) #generating a mask to apply NAs to the final map layer. 	
	Tn <- newdata  #directly allocated all the cells at their most favourable locations according to TTP
	#Tnt <- Tn
	########################################
	runname <- paste0(mtype,"_lu_", (startyear + i),"_run_",loopval,"_of_",numruns,"_m",m,"_r",7-j)
	filepng <- paste0("OUTPUTS/",runname,"_terra.png")
	#filenasc <- paste("lu", (startyear + i), ".asc", sep="")
	
	png(filename = filepng,width = 1200, height = 1200, bg="white") #Plot each land use map to be able to make an animation.
	#plot(T1)
	#plot(T1,breaks=breakpoints,col=colors,main=(startyear+ii))
	plot(Tn,main=(startyear+i))
	dev.off()
		
	#removeTmpFiles(h=5)
	Sys.sleep(0.1)
	setTxtProgressBar(pb, i)
	print (paste("FINISHED: landuse simulation for ", (startyear + i), sep=""))
	writeRaster(Tn,paste0("OUTPUTS/",runname,"_",startyear+i,".tif"),overwrite=T)
	print(proc.time() - ptm)
	#Tn <- raster(Tn)
}
close(pb)
#writeRaster(Tn,paste0("OUTPUTS/",runname,"_raster.tif"),overwrite=T)
writeRaster(Tn,paste0("OUTPUTS/",runname,".tif"),overwrite=T)


	
#######BLOCK TO CALCULATE MAP STATS ON THE FLY
				
if (statsvals > 0) {
	print("creating map statistics table and running map statistics calculation...")
	statline <- mapcompare(Tn,r1,runname)
	statlist <- rbind(statlist,statline)   #add to the statslist file, initially empty, then updated as each map is produced
	} else {
	print("map statistics not required...")
}
		
#######END BLOCK TO CALCULATE MAP STATS ON THE FLY
sim <- Tn  
writeRaster(sim,paste0("OUTPUTS/",mtype,"_sim",endyear,"_run_",loopval,"_of_",numruns,".tif"),overwrite=T)
png(filename = paste("OUTPUTS/lu_",startyear,".png",sep=""),width = 1200, height = 1200, bg="white") #Plot initial map also to animations file
plot(lu1, main = paste("initial land use map - ",startyear,sep=""))  #plot the start map  ## RH MODIFICATION (correcting old error by RH)
dev.off()

sss <- c(lu1,sim)
#sss <- stack(lu1,sim)

names(sss) <- c(paste0("landuse_", startyear),paste0("landuse_", endyear,"_",mtype))
#name1 <- paste("initial land use ",startyear,sep="")
#name1 <- paste("land use in ",startyear,sep="")
#name2 <- paste("simulated land use ",endyear,sep="")
#names(sss) <- c(name1,name2) 
print(plot(sss)) #without the print command, the output will not appear on the screen
#dev.off()
print (paste("simulation complete at ", Sys.time(), sep=""))
#print(plot(sss,breaks=breakpoints,col=colors))   #same as above, but just in case you have colors and breakpoints defined
#save.image(paste0("OUTPUTS/sim_",loopval,".RData"))  #in case you want the workspace objects from each simulation
}

####whichever route we've taken to get here, we need to output map stats if they are needed
if (statsvals > 0) {
#colnames(statlist) <- c("mapname","hits","misses","false alarms","recall (H/(H+M)","Edge_Density","Edge_DensityNA","Fracdim","Fracdim_NA","TCA","TCA_NA","CL","CL_NA","PD","PD_NA")
colnames(statlist) <- c("mapname","hits","misses","false alarms","recall (H/(H+M)","Edge_Density","Fracdim","TCA","CL","PD")
colnames(refline) <- colnames(statlist)
statlist <- rbind(statlist,refline)
write.table(statlist,paste0("OUTPUTS/statlist_",mtype,".csv"),append = FALSE, col.names = NA,sep=",")
}



# function to get variables needed for simulation
setupSims = function(sample,chrom1,chrom2,pos1,pos2)
	{
	# number of rearrangements per sample
	nEvents = table(sample)
	numEvents = mean(nEvents)
	# gamma density
	library(MASS)
	distParams = fitdistr(nEvents,"gamma")
	# probability that rearrangement is on same chromosome
	probSameChrom = mean(chrom1==chrom2)
	# mean distances between breakpoints if on same chromosome
	sameChrom = which(chrom1==chrom2)
	pos = sapply(unique(chrom1[sameChrom]),FUN=function(x) sameChrom[which(chrom1[sameChrom]==x)],simplify=FALSE)
	meanDists = sapply(pos,FUN=function(x) mean(abs(pos2[x]-pos1[x])))
	# chromosome lengths
	lims = getChromInfo()
	chromLengths = lims[2,]-lims[1,]
	return(list(numEvents=numEvents,
		nSamples=length(unique(sample)),
		distrParams=distParams,
		probSameChrom=probSameChrom,
		meanDists=meanDists,
		lims=lims,
		chromLengths=chromLengths))
	}

# funcion to simulate a breakpoint
simulateBreakpoint = function(probSameChrom,chromLengths,meanDists,lims)
	{
	# is it on same chromsome?
	sameChrom = rbinom(1,1,probSameChrom)
	nDraw = ifelse(sameChrom==1,1,2)
	# draw chromosomes
	chroms = sample(colnames(lims),nDraw,prob=chromLengths)
	# get positions
	if(length(chroms)==1)
		{
		SNVlength = round(rexp(1,rate=1/meanDists[chroms]))
		while(SNVlength>lims[2,chroms])
			{
			SNVlength = round(rexp(1,rate=1/meanDists[chroms]))
			}
		start = round(runif(1,lims[1,chroms],lims[2,chroms]-SNVlength))
		pos = c(start,start+SNVlength)
		chroms = c(chroms,chroms)
		} else {
		pos = sapply(chroms,FUN=function(x) round(runif(1,lims[1,x],lims[2,x])))
		}
	out = c(chroms,pos)
	names(out) = c("chrom1","chrom2","pos1","pos2")
	out
	}

# function to simulate a single sample
singleSample = function(probSameChrom,chromLengths,meanDists,lims,distParams)
	{
	#nBreaks = round(runif(1,range(nEvents)[1],range(nEvents)[2])) # uniform
	#nBreaks = max(1,round(rexp(1,1/numEvents))) # exponential
	nBreaks = max(1,round(rgamma(1,shape=distParams$estimate["shape"],rate=distParams$estimate["rate"]))) # gamma
	#print(nBreaks)
	res = replicate(nBreaks,simulateBreakpoint(probSameChrom,chromLengths,meanDists,lims))
	res = rbind(runif(1,0,10000000),res)
	rownames(res)[1] = "Sample"
	res
	}

# function to simulate a single cohort
singleCohort = function(nSamps,probSameChrom,chromLengths,meanDists,lims,distParams,windowsGR,windowsNames)
	{
	res = replicate(nSamps,singleSample(probSameChrom,chromLengths,meanDists,lims,distParams))
	SUBRES <<- res
	res = do.call(cbind,res)
	# GR for result
	samples = rep(res["Sample",],2)
	GRstring = paste0(c(res["chrom1",],res["chrom2",]),":",
			c(res["pos1",],res["pos2",]),"-",
			c(res["pos1",],res["pos2",])) 
	checkGR = grepl("[[:alnum:]]+[:]{1}[[:digit:]]+-{1}[[:digit:]]+$",GRstring)
	if(!all(checkGR)) ERRORGR <<- checkGR
	GR = as(GRstring,"GRanges")
	# overlaps between windows and simulation
	overlaps = findOverlaps(windowsGR,GR)
	# counts of samples  and rearrangements
	counts = table(overlaps@from,samples[overlaps@to])
	sampleCounts = apply(counts,MARGIN=1,FUN=function(x) sum(x!=0))
	# fill in missing
	fullCounts = rep(0,length(windowsNames))
	fullCounts[as.integer(names(sampleCounts))] = sampleCounts
	names(fullCounts) = windowsNames
	if(!is.numeric(fullCounts)) SAVECOHORT <<- list(res=res,samples=samples,GR=GR,overlaps=overlaps,counts=counts,sampleCounts=sampleCounts,fullCounts=fullCounts)
	fullCounts
	}

# function to run a simulation
runSim=function(nCohort,
		sample,chrom1,chrom2,pos1,pos2,
		windowSize=1000000,windowGap=100000,method="slide",
		outDir=getwd(),pmeth="chromosome")
	{
	# setup needed variables
	setup = setupSims(sample,chrom1,chrom2,pos1,pos2)
	SETUP<<-setup
	# windows
	windows = getWindows(windowSize=windowSize,
			windowGap=windowGap,
			method=method)
	WINDOWS<<-windows
	# get simulation results
	res = replicate(nCohort,singleCohort(nSamps=setup$nSamples,
		probSameChrom=setup$probSameChrom,
		chromLengths=setup$chromLengths,
		meanDists=setup$meanDists,
		lims=setup$lims,
		distParams=setup$distrParams,
		windowsGR=windows$windowsGR,
		windowsNames=windows$windows))
	RES<<-res
	resChroms = sapply(rownames(res),FUN=function(x) strsplit(x,split=":")[[1]][1])
	if(pmeth=="chromosome")
		{
		# get P values for each chromosome
		Ps = sapply(1:max(res), FUN=function(x)
			{
			sapply(unique(resChroms), FUN=function(y)
				{
				index = which(resChroms==y)
				(sum(res[index,]>=x)+1)/(prod(dim(res[index,]))+1)
				})
			})
		colnames(Ps) = 1:max(res)
		write.csv(Ps,file=paste0(outDir,"/simulatedP-",method,"-",nCohort,"-chrom.csv"))
		sampleSizes=nCohort*sapply(unique(resChroms),FUN=function(x) length(which(resChroms==x)))
		} else {
		# get P values for each window
		#Ps = apply(res,MARGIN=1,FUN=function(x)
		#	{
		#	sapply(1:max(res), FUN=function(y) (sum(x>=y)+1)/(length(x)+1))
		#	})
		#names(Ps) = 1:max(res)
		Ps = apply(res,MARGIN=1,FUN=function(x)
			{
			tmp = table(x)
			tmp = (rev(cumsum(rev(tmp)))+1)/(sum(tmp)+1)
			})		
		save(list="Ps",file=paste0(outDir,"/simulatedP-",method,"-",nCohort,"-window.Rdata"))
		sampleSizes = rep(nCohort,length(windows))
		names(sampleSizes) = colnames(Ps)
		}
	return(list(Ps=Ps,n=sampleSizes))
	}

#Error in f(init, x[[i]]) : non-numeric argument to binary operator
#In addition: Warning messages:
#1: In densfun(x, parm[1], parm[2], ...) : NaNs produced
#2: In mclapply(cohorts, FUN = function(x) { :
#  scheduled cores 2, 6 encountered errors in user code, all values of the jobs will be affected


#> RES[[2]]
#[1] "Error in asMethod(object) : \n  The character vector to convert to a GRanges object must contain\n  strings of the form \"chr1:2501-2800\" or \"chr1:2501-2800:+\" (\"..\" being\n  also supported as a separator between the start and end positions).\n  Strand can be \"+\", \"-\", \"*\", or missing.\n"
#attr(,"class")
#[1] "try-error"
#attr(,"condition")
#<simpleError in asMethod(object): The character vector to convert to a GRanges object must contain
#  strings of the form "chr1:2501-2800" or "chr1:2501-2800:+" (".." being
#  also supported as a separator between the start and end positions).
#  Strand can be "+", "-", "*", or missing.>

# function to run a simulation
runSimWithTest=function(nCohort,testVals,
		sample,chrom1,chrom2,pos1,pos2,
		windowSize=1000000,windowGap=100000,method="slide",
		outDir=getwd(),pmeth="chromosome",doParallel=FALSE)
	{

	# setup needed variables
	setup = setupSims(sample,chrom1,chrom2,pos1,pos2)
	SETUP<<-setup
	# windows
	windows = getWindows(windowSize=windowSize,
			windowGap=windowGap,
			method=method)
	WINDOWS<<-windows
	if(length(testVals)!=length(windows$windows)) stop("Wrong testVals size")
	# get simulation results
	if(!doParallel)
		{
		# non-parallel
		res = matrix(0,ncol=length(windows$windows),nrow=2)
		for(i in 1:nCohort)
			{
			print(i)
			# run sim
			thisRes = singleCohort(nSamps=setup$nSamples,
				probSameChrom=setup$probSameChrom,
				chromLengths=setup$chromLengths,
				meanDists=setup$meanDists,
				lims=setup$lims,
				distParams=setup$distrParams,
				windowsGR=windows$windowsGR,
				windowsNames=windows$windows)
			if(!is.numeric(thiRes)) ERRORRES <<- thisRes
			# add to succeses
			res[1,] = res[1,]+(thisRes>=testVals)
			# add to total
			res[2,] = res[2,]+1 
			}
		} else {
		# parallel
		# setup up parallel params
		require(parallel)
		ncores = detectCores()
		binSize = ceiling(nCohort/ncores)
		cohorts = 1:nCohort
		cohorts = split(cohorts, ceiling(seq_along(cohorts)/binSize))
		# get simulation results
		res = mclapply(cohorts,FUN=function(x)
			{
			outRes = matrix(0,ncol=length(windows$windows),nrow=2)
			for(i in x)
				{
				print(i)
				# run sim
				thisRes = singleCohort(nSamps=setup$nSamples,
					probSameChrom=setup$probSameChrom,
					chromLengths=setup$chromLengths,
					meanDists=setup$meanDists,
					lims=setup$lims,
					distParams=setup$distrParams,
					windowsGR=windows$windowsGR,
					windowsNames=windows$windows)
				# add to successes
				outRes[1,] = outRes[1,]+(thisRes>=testVals)
				# add to total
				outRes[2,] = outRes[2,]+1 
				}
			outRes
			},mc.cores=ncores)
		RES <<- res
		save(list="res",file=paste0(outDir,"tmp-simRes-",nCohort,"-",method,".Rdata"))
		# combine parallel results
		add <- function(x) Reduce("+", x)
		res=add(res)
		}
	res = res+1 # k+1/n+1 for monte carlo p values
	Ps = res[1,]/res[2,]
	return(Ps=Ps)
	}

# convert a list of counts to p values for each count value
countsToPs = function(counts, # counts in each window (vector)
			MCsims) # MC simulations of each window (list)
	{
	sapply(1:length(counts), FUN=function(x)
		{
		if(counts[x]==0) return(1)
		if(counts[x]>max(as.numeric(names(MCsims[[x]])))) return(min(MCsims[[x]]))
		MCsims[[x]][paste0(counts[x])]
		})
	}



# function to convert monte carlo Ps to Qs
# see Sandve et al, 2011. Sequential Monte Carlo multiple testing
Qval = function(Ps) # p values for each test (each window)
	{
	Ps = Ps[which(!is.na(Ps))]
	n = length(Ps)
	# null proportion estimate
	nullProp = min(1,(2/n)*sum(Ps))
	# ordered p values
	index = 1:length(Ps)
	newOrder = order(Ps)
	ordP = Ps[newOrder]
	index = index[newOrder]	
	# q estimate
	qs = sapply(1:length(ordP),FUN=function(x)
		{
		min(n*nullProp*(ordP[x:n]/(x:n)))
		})
	# reorder qs
	qs = qs[order(index)]
	return(qs)
	}




# function to count rearrangements in genomic windows
countRearrangements = function(sample, # vector of samples for each rearrangement
                               chrom1, # vector of 1st partner chromosomes
                               start1, # vector of 1st partner start positions
                               end1, # vector of 1st partner end positions (can be same as start)
                               chrom2, # vector of 2nd partner chromosomes
                               start2, # vector of 2nd partner start positions
                               end2, # vector of 2nd partner end positions (can be same as start)
                               windowsGR, # GRanges for chromosome windows
                               windowsNames # windows names
                               )
  {
  print("count rearrangements")
  # BRASS genomic ranges
  samples = rep(sample,2)
  brassGranges = as(paste0(c(chrom1,chrom2),":",
                           c(start1,start2),"-",
                           c(end1,end2)),"GRanges")
  # counts
  overlaps = findOverlaps(windowsGR,brassGranges)
  counts = table(overlaps@from,samples[overlaps@to])
  sampleCounts = apply(counts,MARGIN=1,FUN=function(x) sum(x!=0))
  sampleReference = apply(counts,MARGIN=1,FUN=function(x) paste(colnames(counts)[which(x!=0)],collapse=";"))
  # fill in missing counts
  fullCounts = rep(0,length(windowsNames))
  fullCounts[as.integer(names(sampleCounts))] = sampleCounts
  names(fullCounts) = windowsNames
  # fill in missing references
  fullRefs = rep("",length(windowsNames))
  fullRefs[as.integer(names(sampleReference))] = sampleReference
  names(fullRefs) = windowsNames
  # proportion of translocations
  dataGR1 = as(paste0(chrom1,":",start1,"-",end1),"GRanges")   
  dataGR2 = as(paste0(chrom2,":",start2,"-",end2),"GRanges")  
  matchingProp = sapply(windowsGR,FUN=function(x)
	{
	overlapping = unique(c(findOverlaps(dataGR1,x)@from,findOverlaps(dataGR2,x)@from))
	matching = chrom1[overlapping]!=chrom1[overlapping]
	sum(matching)/length(matching)
	})
  names(matchingProp) = windowsNames
  # tidy
  rownames(counts) = windowsNames[as.numeric(rownames(counts))] 
  return(list(sampleCounts=fullCounts,dataGR=brassGranges,
	sampleReference=fullRefs,detailedCounts=counts,
	translocProp=matchingProp))
}

# save res
saveRes = function(sampleCounts,sampleReference,windows,cutoff=NULL,outDir=getwd())
  {
  print("important regions")
  out = cbind(windows,sampleCounts,sampleReference)
  if(is.null(cutoff)) cutoff = ceiling(length(sampleCounts)*0.07) # more than 7% of samples
  toSave = out[which(as.numeric(out[,2])>cutoff),]
  colnames(toSave)[1] = "regions"

  write.csv(toSave,paste0(outDir,"/recurrentRearrangements.csv"),row.names=FALSE,quote=FALSE)
  return(toSave)
  }


# function to generate windows and count rearrangements only
getRearrCounts = function(sample,chrom1,start1,end1,
                                chrom2,start2,end2,
                                cutoff=5,outDir=getwd(),
                                windowSize=1000000,windowGap=100000,windowMethod="slide"
                               )
	{
	# get genomic windows
	windows = RETREAD:::getWindows(windowSize=windowSize,windowGap=windowGap,method=windowMethod)
	# count rearrangments in windows
	rearrangementCounts = RETREAD:::countRearrangements(sample=sample, # vector of samples for each rearrangement
                chrom1=chrom1, # vector of 1st partner chromosomes
                start1=start1, # vector of 1st partner start positions
                end1=end1, # vector of 1st partner end positions (can be same as start)
                chrom2=chrom2, # vector of 2nd partner chromosomes
                start2=start2, # vector of 2nd partner start positions
                end2=end2, # vector of 2nd partner end positions (can be same as start)
                windowsGR=windows$windowsGR, # GRanges for chromosome windows
                windowsNames=windows$windows # windows names
		)
	# save interesting windows
	toPlot = RETREAD:::saveRes(sampleCounts=rearrangementCounts$sampleCounts,
		sampleReference=rearrangementCounts$sampleReference,
                   windows=windows$windows,
                   cutoff=cutoff,
                   outDir=outDir)
	return(list(windows=windows,rearrCounts=rearrangementCounts,info=toPlot))
	}

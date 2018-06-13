# tiling windows
getWindowsTile = function(lims,windowSize)
  {
  print("tiling windows")
  windowsTile = unlist(sapply(colnames(lims),FUN=function(x)
    {
    splits = as.integer(seq(from=lims[1,x],to=lims[2,x],by=windowSize))
    starts = splits[-length(splits)]
    ends = splits[-1]
    paste0(x,":",starts,"-",ends)
    }))
  return(windowsTile)
  }

# sliding windows
getWindowsSlide = function(lims,windowGap,windowSize)
  {
  print("sliding windows")
  windowsSlide = unlist(sapply(colnames(lims),FUN=function(x)
    {
    starts = as.integer(seq(from=lims[1,x],to=lims[2,x]-windowSize,by=windowGap))
    ends = starts+windowSize
    paste0(x,":",starts,"-",ends)
    }))
  return(windowsSlide)
  }

# function to get genomic windows
getWindows = function(windowSize=1000000,windowGap=100000,method="slide")
  {
  print("get windows")
  lims = getChromInfo()
  if(method=="slide")
    {
    windows = getWindowsSlide(lims,windowGap,windowSize)
    } else {
    windows = getWindowsTile(lims,windowSize)
    }
  windowsGranges = as(windows,"GRanges")
  return(list(windows=windows,windowsGR=windowsGranges,lims=lims))
  }

# count number of chromothriptic regions overlapping genomic windows
countChromothripsis = function(chromothripticRes,windowsGR,windowsNames)
  {
  print("count chromothripsis")
  chromothriptic = do.call(rbind,sapply(chromothripticRes,FUN=function(x) do.call(rbind,x)))
  chromothripGR = as(paste0(chromothriptic[,"chrom"],":",
                            as.integer(chromothriptic[,"windowStarts"]),"-",
                            as.integer(chromothriptic[,"windowEnds"])),
                     "GRanges")
  # chromothripsis counts in windows
  overlaps = findOverlaps(windowsGR,chromothripGR)
  counts = table(overlaps@from,chromothriptic[,"samp"][overlaps@to])
  chromoSampleCounts = apply(counts,MARGIN=1,FUN=function(x) sum(x!=0))
  chromoFullCounts = rep(0,length(windowsGR))
  
  chromoFullCounts[as.integer(names(chromoSampleCounts))] = chromoSampleCounts

  # expected number of chromothripses per chromosome
  combined = unique(chromothriptic[,1:2])
  # check if chromothripsis is randomly distributed across chromosomes
  chromoCounts = sapply(c(1:22,"X","Y"),FUN=function(x) sum(combined[,2]==x))
  chromInfo = getChromInfo()
  chromLengths = sapply(c(1:22,"X","Y"),FUN=function(x) abs(diff(chromInfo[,x])))
  # goodness of fit
  test = chisq.test(chromoCounts,p=chromLengths/sum(as.numeric(chromLengths))) # not significantly different
  expected = sum(chromoCounts)*(chromLengths/sum(as.numeric(chromLengths)))
  # tidy
  rownames(counts) = windowsNames[as.numeric(rownames(counts))] 
  return(list(chromoCounts=chromoFullCounts,expected=expected,detailedCounts=counts))
}

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

# function to plot just sample counts
plotSampleCounts=function(counts, # rearrangment counts
                          colourGroup, # one per interesting window, colour group
                          groupLabels, # one per group, labels
                          groupColours, # one per group, colours
                          important, # interesting results (toSave)
                          cutoff=cutoff
                          )
  {
  print("plot rearrangement counts")
  # plot sample counts
  par(mar=c(0,4,0,2))
  plot(counts,type="l",xaxt="n",ylab="# Rearranged",bty="n",col=
"gray37")
  abline(h=cutoff,lty=2)
  for(i in unique(colourGroup))
    {
    index = which(colourGroup==i)
    regions = important[index,1]
    index = which(names(counts)%in%regions)
    x = min(index)+max(sapply(strsplit(groupLabels[i],split="\n")[[1]],FUN=nchar))*18
    y = max(counts[index])+(diff(range(counts))*0.025)
    text(x=x,y=y,labels=groupLabels[i],col=groupColours[i])
    }
  }

# function to plot jst chromothriptic counts
plotChromothripsisCounts=function(counts, # chromothriptic counts
                                expected,  # expected number per chromosome
				windowNames,
				plotExpectedChromothriptic=TRUE
				)
  {
SAVECOUNTS <<- counts
SAVEEXPECTED <<- expected
  print("plot chromothriptic counts")
  # plot chromothriptic counts
  par(mar=c(0,4,0,2))
  plot(counts,type="l",lwd=2,col="gray37",xaxt="n",ylab="# Chromothriptic",bty="n")
  if(plotExpectedChromothriptic)
	{
  	chroms = sapply(windowNames,FUN=function(x) strsplit(x,split="[:]|-")[[1]][1])
  	for(i in 1:length(expected))
                {
		index = which(chroms==names(expected)[i])
                end = max(index)
                start=min(index)
                lims = c(start,end)
                lines(x=lims,y=rep(expected[i],2),lty=2)
                }
	} 
 }

# function to plot just chromosome bar
plotChromsomeBar=function(windows)
  {
  print("plot chromosomes bar")
  # plot chromosome bar
  par(mar=c(0,4,0,2))
  rangeExtra = 0.4
  textSep = 0.25
  plot(NA,xlim=c(1,length(windows)),ylim=c(0-rangeExtra,1+rangeExtra),bty="n",yaxt="n",xaxt="n",ylab=NA)
  chroms = sapply(windows,FUN=function(x) strsplit(x,split="[:]")[[1]][1])
  colours = c("white","gray")
  plotCol = colours[1]
  for(i in unique(chroms))
    {
    plotCol=colours[which(colours!=plotCol)]
    lims = range(which(chroms==i))
    polygon(x=c(lims,rev(lims)),y=c(0,0,1,1),col=plotCol,border="gray56")
    text(x=mean(lims),y=ifelse(plotCol=="white",0-textSep,1+textSep),labels=i,cex=0.8)
    }
  }

# function to plot a single rearrangment
getCurve = function(start,end,colour="black",squishFactor=1)
  {
  if(start==end)
    {
    lines(c(start,end),c(0,1),col=colour,lwd=2)
    return()
    }
  vals = sort(c(start,end))
  x = seq(from=vals[1],to=vals[2],length.out=100)
  xNorm = (median(x)-x)/abs(diff(range(x)))
  y = xNorm^2
  y = y/max(y)
  y = y*squishFactor
  y = y+(1-max(y))
  print(c(range(x),range(y),colour))
  lines(x,y,type="l",col=colour,lwd=2)
  }

# function to plot just rearrangements
plotRearrangements=function(toPlot, # string of windows within which to plot rearrangments (from toSave)
                            dataGR, # data genomic ranges
                            windowsGR, # windows genomic ranges
                            colourGroup, # colour group
                            colours, # colours
                            chrom1, # data chrom1
                            start1, # data start1
                            end1, # data end1 
                            chrom2, # data chrom2
                            start2, # data start2
                            end2 # data end2
                            )
  {
  print("plot rearrangements")
  par(mar=c(0,4,0,2))
  plot(NA,xlim=c(1,length(windowsGR)),ylim=0:1,bty="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
  # rearrangements to plot
  interesting = as(toPlot[,1],"GRanges")
  overlaps1 = findOverlaps(interesting,dataGR[1:(length(dataGR)/2)])
  overlaps2 = findOverlaps(interesting,dataGR[((length(dataGR)/2)+1):length(dataGR)])
  info = cbind(c(overlaps1@to,overlaps2@to),
               c(colours[colourGroup[overlaps1@from]],colours[colourGroup[overlaps2@from]]),
               c(colourGroup[overlaps1@from],colourGroup[overlaps2@from]))
  overlaps = unique(info)
  overlaps = overlaps[rev(order(overlaps[,3])),]
  windowJoins = apply(overlaps,MARGIN=1,FUN=function(x)
  {
    print(x)
    i = as.integer(x[1])
    window1 = findOverlaps(windowsGR,as(paste0(chrom1[i],":",start1[i],"-",end1[i]),"GRanges"))@from
    window2 = findOverlaps(windowsGR,as(paste0(chrom2[i],":",start2[i],"-",end2[i]),"GRanges"))@from
    if(length(window1)==0|length(window2)==0)
    {
      print(paste0("Missing: ",x))
      return(c(window1,window2,NA))
    }
    squishFactor = (abs(diff(c(window1,window2)))+100)/length(windowsGR)
    if(squishFactor<0.01) squishFactor=0.8
    print(c(window1,window2,x[2],squishFactor))
    getCurve(window1,window2,col=x[2],squishFactor=squishFactor)
    return(c(window1,window2))
  })
  }


# function to plot everything together
plotCounts = function(rearrangementCounts, # rearrangment counts
                      colourGroup, # one per interesting window, colour group
                      groupLabels, # one per group, labels
                      toPlot, # which windows to plot (above cutoff)
                      chromothripticCounts, # chromothripsis counts
                      windows, # window labels e.g. X:1234-5678
                      dataGR, # data genomic ranges
                      windowsGR, # windows genomic ranges
                      chrom1=chrom1, # data chrom1
                      start1=start1, # data start1
                      end1=end1, # data end1 
                      chrom2=chrom2, # data chrom2
                      start2=start2, # data start2
                      end2=end2, # data end2
                      cutoff=cutoff,
		      lims = NULL,
		      method="slide",
			plotExpectedChromothriptic=TRUE,
                      groupColours =  c("lightblue","dodgerblue3",
                                        "lightgreen","green4",
                                        "lightpink","red2",
                                        "peachpuff2","orange",
                                        "plum3","purple3") # one per group, indicating colour
                      )
  {
  print("begin plotting")
  layout(matrix(c(1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  1,1,1,
                  2,2,2,
                  2,2,2,
                  2,2,2,
                  3,3,3,
                  4,4,4,
                  4,4,4,
                  4,4,4
  ),byrow=TRUE))
  if(method=="slide")
	{
	# sliding window
	plotSliding(rearrangeCounts=rearrangementCounts,
		chromoCounts=chromothripticCounts$chromoCounts,
		chromoExpected = chromothripticCounts$expected,
		chromEnds=getPosVec(lims),
		toPlot=toPlot,
		dataGR=dataGR,
		colourGroup=colourGroup,
		groupLabels=groupLabels,
		groupColours=groupColours,
		chrom1=chrom1,
		chrom2=chrom2,
		pos1=start1,
		pos2=start2,
		cutoff=cutoff)
	} else {
	# tiling window
  	plotSampleCounts(counts=rearrangementCounts, # rearrangment counts
                   colourGroup=colourGroup, # one per interesting window, colour group
                   groupLabels=groupLabels, # one per group, labels
                   groupColours=groupColours, # one per group, colours
                   important=toPlot,
                   cutoff=cutoff)
  	plotChromothripsisCounts(counts=chromothripticCounts$chromoCounts,expected=chromothripticCounts$expected,windowNames=windows,plotExpectedChromothriptic=plotExpectedChromothriptic)
	plotChromsomeBar(windows=windows)
	plotRearrangements(toPlot=toPlot, # string of windows within which to plot rearrangments (from toSave)
                     dataGR=dataGR, # data genomic ranges
                     windowsGR=windowsGR, # windows genomic ranges
                     colourGroup=colourGroup, # colour group
                     colours=groupColours, # colours
                     chrom1=chrom1, # data chrom1
                     start1=start1, # data start1
                     end1=end1, # data end1 
                     chrom2=chrom2, # data chrom2
                     start2=start2, # data start2
                     end2=end2)
	}
}

# function to run everything
recRearrangePipeline = function(sample,chrom1,start1,end1,
                                chrom2,start2,end2,
                                chromothripticRes,
                                cutoff=5,outDir=getwd(),
                                windowSize=1000000,windowGap=100000,windowMethod="slide",
                                colourGroup=NULL,groupLabels=NULL,
                                cosmicGenes=NULL,
				plotExpectedChromothriptic=TRUE,
                                groupColours=c("lightblue","dodgerblue3",
                                               "lightgreen","green4",
                                               "lightpink","red2",
                                               "peachpuff2","orange",
                                               "plum3","purple3") # one per group, indicating colour
                                )
  {
  print("run analysis")
  # get genomic windows
  windows = getWindows(windowSize=windowSize,windowGap=windowGap,method=windowMethod)
  # count rearrangments in windows
  rearrangementCounts = countRearrangements(sample=sample, # vector of samples for each rearrangement
                                   chrom1=chrom1, # vector of 1st partner chromosomes
                                   start1=start1, # vector of 1st partner start positions
                                   end1=end1, # vector of 1st partner end positions (can be same as start)
                                   chrom2=chrom2, # vector of 2nd partner chromosomes
                                   start2=start2, # vector of 2nd partner start positions
                                   end2=end2, # vector of 2nd partner end positions (can be same as start)
                                   windowsGR=windows$windowsGR, # GRanges for chromosome windows
                                   windowsNames=windows$windows # windows names
    )
  # count chromothipses in windows
  chromoCounts = countChromothripsis(chromothripticRes=chromothripticRes,
                                     windowsGR=windows$windowsGR,
				     windowsNames=windows$windows)
  # save interesting windows
  toPlot = saveRes(sampleCounts=rearrangementCounts$sampleCounts,
		sampleReference=rearrangementCounts$sampleReference,
                   windows=windows$windows,
                   cutoff=cutoff,
                   outDir=outDir)
  TOPLOT <<- toPlot
  WINDOWS <<- windows
  SAMPLECOUNTS <<- rearrangementCounts
  # setup groups
  if(is.null(colourGroup))
    {
    #colourGroup = 1:nrow(toPlot)
    #groupLabels = paste0(sapply(toPlot[,1],FUN=function(x) strsplit(x,split="[:]")[[1]][1]),1:nrow(toPlot))
    tmp = autoLabel(toPlot,cosmicGenes)
    colourGroup = tmp$groups
    groupLabels = tmp$labels
    }
  if(length(groupLabels)>length(groupColours))
    {
    groupColours = c(groupColours,
                     rainbow(length(groupLabels)-length(groupColours)))
    }
  COLOURS <<- groupColours
  GROUPS <<- colourGroup
  LABELS <<- groupLabels
  # plot everything
  plotCounts(rearrangementCounts=rearrangementCounts$sampleCounts, # rearrangment counts
              colourGroup=colourGroup, # one per interesting window, colour group
              groupLabels=groupLabels, # one per group, labels
              toPlot=toPlot, # which windows to plot (above cutoff)
              chromothripticCounts=chromoCounts, # chromothripsis counts
              windows=windows$windows, # window labels e.g. X:1234-5678
              dataGR=rearrangementCounts$dataGR, # data genomic ranges
              windowsGR=windows$windowsGR, # windows genomic ranges
              chrom1=chrom1, # data chrom1
              start1=start1, # data start1
              end1=end1, # data end1 
              chrom2=chrom2, # data chrom2
              start2=start2, # data start2
              end2=end2, # data end2
              cutoff=cutoff, # cutoff
	      lims=windows$lims,
	      method=windowMethod,
              groupColours =  groupColours, # one per group, indicating colour
		plotExpectedChromothriptic=plotExpectedChromothriptic
      )
  toPlot = cbind(toPlot,unlist(sapply(unique(colourGroup),FUN=function(x) rep(groupLabels[x],length(which(colourGroup==x))))))
  colnames(toPlot)[ncol(toPlot)] = "group"
  return(list(windows=windows,
              rearrangeCounts=rearrangementCounts,
              chromoCounts=chromoCounts,
              toSave=toPlot,
              group = colourGroup,
              labels = groupLabels
              ))
  }

# function to automatically label regions
autoLabel = function(toSave,cosmicGenes,doChromBand=FALSE)
{
  # get overlapping groups
  GR = as(toSave[,1],"GRanges")
  # combine regions
  combReg = combineRegions(regions=t(sapply(paste0(toSave[,1]),FUN=function(x) strsplit(x,split="[:]|-")[[1]])),
		chromCol=1,startCol=2,endCol=3)
  combinedGR = as(paste0(combReg[,1],":",combReg[,2],"-",combReg[,3]),"GRanges")
  # groups from combined regions
  overlaps = findOverlaps(GR,combinedGR)
  groups = overlaps@to
  # annotate regions with genes of interest
  library(biomaRt)
  genes = sapply(toSave[,1],FUN=function(x) 
    {
    info = strsplit(x,split=c("-|[:]"))[[1]]
    getOverlapGenes(info[1],info[2],info[3])
    })
  # group labels per region
  groupLabels = sapply(genes["geneName",],FUN=function(x) 
    {
    index = which(cosmicGenes%in%x)
    if(length(index)>0) return(paste(cosmicGenes[index],collapse=";"))
    return(NA)
    })
  # group labels per group
  groupLabels = sapply(unique(groups),FUN=function(x)
    {
    match = groupLabels[which(groups==x)]
    naIndex = which(is.na(match))
    if(length(naIndex)==0)
      {
      out = paste(unique(match),collapse=";")
      return(paste(unique(strsplit(out,split=";")[[1]]),collapse=";"))
      }
    if(length(naIndex)!=length(match))
      {
      out=paste(unique(match[-naIndex]),collapse=";")
      return(paste(unique(strsplit(out,split=";")[[1]]),collapse=";"))
      }
    NA
    })
 # fill in missing
  index = which(is.na(groupLabels))
  if(length(index)>0)
	{
	groupLabels[index] = sapply(index,FUN=function(x)
		{
		info = sapply(toSave[which(groups==x),1],FUN=function(y) strsplit(y,split="[:]|-")[[1]],simplify=FALSE)
		BLAH <<- info
		info = do.call(rbind,info)
		chrom=info[1,1]
		start=min(info[,2:3])
		end=max(info[,2:3])
		start = round(as.numeric(start)/1000000,1)
		end = round(as.numeric(end)/1000000,1)
		if(doChromBand)
			{
			chromBand = getChromBand(chrom=chrom,start=start,end=end)
			paste0(chromBand,"\n[",chrom,":",start,"-",end,"Mb]")
			} else {
			paste0("[",chrom,":",start,"-",end,"Mb]")
			}
		})
	}
  return(list(groups=groups,labels=groupLabels))
}

# function to get chromsome band
getChromBand = function(chrom,start,end)
	{
	library(biomaRt)
	ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
	bandInfo = getBM(attributes=c("band"),
		filters=c("chromosome_name","start","end"),
		values=list(chromosome_name=chrom,
			start=start,
			end=end),
		mart=ensembl)
	return(paste0(chrom,bandInfo$band[1]))
	}


# function to get genes that overlap a region
getOverlapGenes = function(chr,start,end)
{
  # get gRanges from this DMR
  ranges = gsub(" ","",paste0(chr,":",start,":",end))
  library(dplyr)
  gRanges = sapply(ranges, function (x) {res=strsplit(x, ':')}) %>%
    unlist %>%
    as.numeric %>%
    matrix(ncol=3, byrow=T) %>%
    as.data.frame %>%
    dplyr::select(chrom=V1, start=V2, end=V3) %>%
    mutate(chrom=paste0('chr', chrom)) %>%
    makeGRangesFromDataFrame
  # get Hsapiens genes
  library(Homo.sapiens)
  genesRanges = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  # overlaps between genes and DMR
  overlap = subsetByOverlaps(genesRanges,gRanges); overlap
  geneStart = overlap@ranges@start
  geneEnd = overlap@ranges@start+overlap@ranges@width
  # get gene names
  library(org.Hs.eg.db)
  print(str(org.Hs.egSYMBOL))
  unmapped = org.Hs.eg.db::org.Hs.egSYMBOL
  mapped = mappedkeys(unmapped)
  genes = unlist(as.list(unmapped[mapped]))
  DMRgenes = genes[which(names(genes)%in%overlap$gene_id)]
  # return
  return(list(geneStart=geneStart,geneEnd=geneEnd,geneName=DMRgenes))
}


# function to combine overlapping regions
combineRegions = function(regions,chromCol=2,startCol=3,endCol=4)
	{
	if(nrow(regions)<2) return(regions)
	doOverlap=TRUE
	while(doOverlap)
		{
		seqStrings = paste0("chr",regions[,chromCol],":",regions[,startCol],"-",regions[,endCol])
		seqGranges = as(seqStrings,"GRanges")
		overlaps = findOverlaps(seqGranges,seqGranges)
		doOverlap = length(overlaps@from)!=nrow(regions)
		if(!doOverlap) break
		overlapInfo = sapply(1:length(seqStrings),FUN=function(x) overlaps@to[which(overlaps@from==x)],simplify=FALSE)
		newRegions = t(sapply(overlapInfo,FUN=function(x) c(
					unique(regions[x,chromCol]),
					min(as.numeric(regions[x,c(startCol,endCol)])),
					max(as.numeric(regions[x,c(startCol,endCol)]))
					)))
		regions = unique(newRegions)
		}
	return(regions)
	}


getPosVec = function(lims)
	{
	chromEnds = cumsum(as.numeric(lims[2,]-lims[1,]))
	names(chromEnds) = colnames(lims)
	return(chromEnds)
	}

# function to plot with sliding windows
plotSliding = function(rearrangeCounts,chromoCounts,chromoExpected,
			chromEnds,toPlot,
			dataGR,colourGroup,groupLabels,groupColours,
			chrom1,chrom2,pos1,pos2,
			rangeExtra = 0.4,textSep = 0.25, cutoff,plotExpectedChromothriptic)
	{
	par(mar=c(0,4,0,2))
	chromLengths = vector(length=length(chromEnds))
	chromLengths[2:length(chromLengths)] = chromEnds[2:length(chromLengths)]-chromEnds[1:(length(chromLengths)-1)]
	chromLengths[1] = chromEnds[1]
	names(chromLengths) = names(chromEnds)
	# convert genomic information to plotting coordinates
	info = sapply(names(rearrangeCounts),FUN=function(x) strsplit(x,split="[:]|-")[[1]])
	info = apply(info,MARGIN=2,FUN=function(x) 
		{
		index = which(names(chromEnds)==x[1])
		chromEnds[index]-(chromLengths[index]-as.integer(x[2:3]))
		})
	info = rbind(info,rearrangeCounts,chromoCounts)
	# plot rearrangement counts
	print("plot rearrangement counts")
	plot(NA,xlim=c(1,max(chromEnds)),ylim=c(0,max(rearrangeCounts)),bty="n",xlab=NA,xaxt="n",ylab="# Rearranged")
	apply(info,MARGIN=2,FUN=function(vec) polygon(x=c(vec[1:2],vec[2:1]),y=c(rep(vec[3],2),0,0),col="gray37",border="gray37"))
	abline(h=cutoff,lty=2)
	for(i in unique(colourGroup))
		{
		index = which(colourGroup==i)
		thisInfo = sapply(toPlot[index,1],FUN=function(x) strsplit(x,split="-|[:]")[[1]],simplify=FALSE)
		thisInfo = do.call(rbind,thisInfo)
		chrom = unique(thisInfo[,1])
		pos = mean(range(as.numeric(thisInfo[,2:3])))
		x = chromEnds[chrom]-(chromLengths[chrom]-pos)
		y = max(as.numeric(toPlot[index,2]))+0.5
		text(x=x,y=y,label=groupLabels[i],col=groupColours[i])
		}
	# plot chromothripsis counts
	print("plot chromothriptic counts")
	plot(NA,xlim=c(1,max(chromEnds)),ylim=c(0,max(chromoCounts)),bty="n",xlab=NA,xaxt="n",ylab="# Chromothriptic")
	apply(info,MARGIN=2,FUN=function(vec) polygon(x=c(vec[1:2],vec[2:1]),y=c(rep(vec[4],2),0,0),col="gray37",border="gray37"))
	if(plotExpectedChromothriptic)
		{
		for(i in 1:length(chromoExpected))
                	{
                	end = chromEnds[i]
                	start=ifelse(i==1,1,chromEnds[i-1])
                	lims = c(start,end)
                	lines(x=lims,y=rep(chromoExpected[i],2),lty=2)
                	}
		}

	# plot chromosomes
	print("plot chromosomes")
	plot(NA,xlim=c(1,max(chromEnds)),ylim=c(0-rangeExtra,1+rangeExtra),bty="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
	colours = c("white","gray")
  	plotCol = colours[1]
	for(i in 1:length(chromEnds))
		{
		plotCol=colours[which(colours!=plotCol)]
		end = chromEnds[i]
		start=ifelse(i==1,1,chromEnds[i-1])
		lims = c(start,end)
		polygon(x=c(lims,rev(lims)),y=c(0,0,1,1),col=plotCol,border="gray56")
    		text(x=mean(lims),y=ifelse(plotCol=="white",0-textSep,1+textSep),labels=names(chromEnds)[i],cex=0.8)
		}
	# plot rearrangements
	print("plot rearrangements")
 	plot(NA,xlim=c(1,max(chromEnds)),ylim=0:1,bty="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
	# rearrangements to plot
  	interesting = as(toPlot[,1],"GRanges")
	overlaps1 = findOverlaps(interesting,dataGR[1:(length(dataGR)/2)])
	overlaps2 = findOverlaps(interesting,dataGR[((length(dataGR)/2)+1):length(dataGR)])
	info = cbind(c(overlaps1@to,overlaps2@to),
	               c(groupColours[colourGroup[overlaps1@from]],groupColours[colourGroup[overlaps2@from]]),
	               c(colourGroup[overlaps1@from],colourGroup[overlaps2@from]))
	overlaps = unique(info)
	overlaps = overlaps[rev(order(overlaps[,3])),]
	# plot the individual rearrangements
	apply(overlaps,MARGIN=1,FUN=function(x)
		{
		index = as.integer(x[1])
		chrom1index = which(names(chromEnds)==chrom1[index])
		chrom2index = which(names(chromEnds)==chrom2[index])
		start = chromEnds[chrom1index]-(chromLengths[chrom1index]-pos1[index])
		end = chromEnds[chrom2index]-(chromLengths[chrom2index]-pos2[index])
		squishFactor = (abs(diff(c(start,end)))+100)/max(chromEnds)
    		if(squishFactor<0.01) squishFactor=0.8
		getCurve(start,end,col=x[2],squishFactor=squishFactor)
		})
	}
	

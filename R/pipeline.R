# function to run analysis and then plot
recRearrangePipeline = function(sample,chrom1,start1,end1,
                                chrom2,start2,end2,
                                cutoff=5,outDir=getwd(),
                                windowSize=1000000,windowGap=100000,windowMethod="slide",
                                colourGroup=NULL,groupLabels=NULL,
                                cosmicGenes=NULL,
                                groupColours=c("lightblue","dodgerblue3",
                                               "lightgreen","green4",
                                               "lightpink","red2",
                                               "peachpuff2","orange",
                                               "plum3","purple3") # one per group, indicating colour
                                )
  {
  print("run analysis")
  # count rearrangements
  data = getRearrCounts(sample=sample,chrom1=chrom1,start1=start1,end1=end1,
                                chrom2=chrom2,start2=start2,end2=end2,
                                cutoff=cutoff,outDir=outDir,
                                windowSize=windowSize,windowGap=windowGap,windowMethod=windowMethod
                               )
  windows = data$windows
  rearrangementCounts = data$rearrCounts
  toPlot = data$info
  TOPLOT <<- toPlot
  WINDOWS <<- windows
  SAMPLECOUNTS <<- rearrangementCounts
  # setup groups
  if(is.null(colourGroup))
    {
    #colourGroup = 1:nrow(toPlot)
    #groupLabels = paste0(sapply(toPlot[,1],FUN=function(x) strsplit(x,split="[:]")[[1]][1]),1:nrow(toPlot))
    tmp = RETREAD:::autoLabel(toPlot,cosmicGenes)
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
  RETREAD:::plotCounts(rearrangementCounts=rearrangementCounts$sampleCounts, # rearrangment counts
              colourGroup=colourGroup, # one per interesting window, colour group
              groupLabels=groupLabels, # one per group, labels
              toPlot=toPlot, # which windows to plot (above cutoff)
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
              groupColours =  groupColours # one per group, indicating colour
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


	

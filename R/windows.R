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

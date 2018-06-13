# load cytoBand file
loadCytoBand = function(file=NULL)
	{
	if(is.null(file))
		{
		tmpEnv = new.env()
		data(list="cytoBand", package='CNomplexity',envir=tmpEnv)
		return(tmpEnv[["cytoBand"]])
		} else {
		return(readFile(file))
		}
	}


# get chromosome lengths
getChromInfo = function(cytoFile=NULL)
	{
	data = loadCytoBand(cytoFile)
	data[,1] = gsub("chr","",data[,1])
	sapply(c(1:22,"X","Y"),FUN=function(x) 
		c(min(data[which(data[,1]==x),2]),
		max(data[which(data[,1]==x),3])))
	}


\name{recRearrangePipeline} 
\alias{recRearrangePipeline} 
\title{RETREAD::recRearrangePipeline} 
\usage{recRearrangePipeline(sample,chrom1,start1,end1,
                                chrom2,start2,end2,
                                cutoff=5,outDir=getwd(),
                                windowSize=1000000,windowGap=100000,windowMethod="tile",
                                colourGroup=NULL,groupLabels=NULL,
                                cosmicGenes=NULL,
                                groupColours=c("lightblue","dodgerblue3",
                                               "lightgreen","green4",
                                               "lightpink","red2",
                                               "peachpuff2","orange",
                                               "plum3","purple3") }
\arguments{
\item{sample}{Vector with sample of real rearrangements, one per rearrangement.}
\item{chrom1}{Vector with chromosome of real rearrangements first partner, one per rearrangement.}
\item{start1}{Vector with start position of real rearrangements first partner, one per rearrangement.}
\item{end1}{Vector with end position of real rearrangements first partner, one per rearrangement.}
\item{chrom2}{Vector with chromosome of real rearrangements second partner, one per rearrangement.}
\item{start2}{Vector with start position of real rearrangements second partner, one per rearrangement.}
\item{end2}{Vector with end position of real rearrangements second partner, one per rearrangement.}
\item{cutoff}{Value above which to plot recurrent rearrangements.}
\item{outDir}{Out directory.}
\item{windowSize}{Window size to split genome into.}
\item{windowGap}{Gap between windows for sliding windows (\code{windowMethod="slide"}).}
\item{windowMethod}{Method to split genome into windows. "slide"=sliding windows, "tile"=tiling windows.}
\item{colourGroup}{Colours for groups. Default to NULL.}
\item{groupLabels}{Labels for groups. Defaults to NULL.}
\item{cosmicGenes}{Genes to highlight on plot.}
\item{groupColours}{Colours for groups. If more groups are identified than specified here, the rest are filled with rainbow.}
}
\value{
List with 5 elements. \code{windows} are the windows that the genome was split into. \code{rearrangeCounts} are the counts of the number of samples with a rearrangment in each window. \code{toSave} is a matrix with details of the windows that have more than \code{cutoff} samples rearranged. \code{group} gives the grouped windows. \code{labels} gives labels for each \code{group}.
} 
\description{
Identify regions that are recurrently rearranged.
}
\details{
Recurrent rearrangement, with no statistical test. Pass \code{rearrangeCounts} to \code{runSimWithTest} to get p-values, and the results of that into \code{Qval} to get q-values.
}


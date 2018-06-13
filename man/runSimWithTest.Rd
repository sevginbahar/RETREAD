\name{runSimWithTest} 
\alias{runSimWithTest} 
\title{RETREAD::runSimWithTest} 
\usage{runSimWithTest(nCohort,testVals,sample,chrom1,chrom2,pos1,pos2,
windowSize=1000000,windowGap=100000,method="slide",
outDir=getwd(),pmeth="chromosome",doParallel=FALSE)}
\arguments{
\item{nCohort}{Number of cohorts to simulate.}
\item{testVals}{Counts of rearrangments in windows from real data set.}
\item{sample}{Vector with sample of real rearrangements, one per rearrangement.}
\item{chrom1}{Vector with chromosome of real rearrangements first partner, one per rearrangement.}
\item{chrom2}{Vector with chromosome of real rearrangements second partner, one per rearrangement.}
\item{pos1}{Vector with position of real rearrangements first partner, one per rearrangement.}
\item{pos2}{Vector with position of real rearrangements second partner, one per rearrangement.}
\item{windowSize}{Size of window used to count rearrangements in real dataset. Defaults to 1000000.}
\item{windowGap}{Gap between windows for sliding windows. Defaults to 100000.}
\item{method}{Sliding ("slide") or tiling ("tile") windows. Default to "tile". Note that the current q-value implementation does not work with sliding windows.}
\item{outDir}{Directory to save to.}
\item{pmeth}{Method to calculate p-values. "chromosome" caluclates p-values per chromosome. "all" calculates p-values for the whole genome. Defaults to "chromosome."}
\item{doParallel}{Whether to perform parallel computation. Defaults to FALSE.}
}
\value{
Vector of p-values.
} 
\description{
Simulates \code{nCohort} new datasets to generate p-values for recurrent rearrangments.
}
\details{
Generates p-values for recurrent rearrangements.
}


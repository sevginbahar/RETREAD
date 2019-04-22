# RETREAD
REcurrenT REArrangement Discovery

For full details of the method and citation, please see [Steele et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1535610819300972).

Dependencies: parallel, GenomicRanges

## Method to identify recurrent rearrangments
Input: bedpe

Simulates n "null" samples of randomly distributed rearrangements, to calculate p-value and q-value for real data. 

getRearrCounts() -> runSimWithTest() -> Qval()

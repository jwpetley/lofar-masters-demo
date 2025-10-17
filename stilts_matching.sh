# Stilts demo
# Documentation here - https://www.star.bris.ac.uk/~mbt/stilts/


topcat -stilts --help # Find all the functions available within stilts

CAT1=./CATALOGUES/LOFAR_EDFN_v1.0.fits
CAT2=./CATALOGUES/master_euclid_edfn-result.fits

CAT3=./CATALOGUES/lofar_euclid_matched.fits


## Plotting functions

# Plot sky positions
# topcat -stilts plot2sky layer1=mark in1=$CAT1 lon1=RA lat1=DEC out=./PLOTS/EDFN_lofar.pdf

# topcat -stilts plot2sky layer1=mark in1=$CAT2 lon1=right_ascension lat1=declination \
# layer2=mark in2=$CAT1 lon2=RA lat2=DEC \
# out=./PLOTS/EDFN_lofar.pdf


# # Plot histogram of fluxes
# topcat -stilts plot2plane layer=histogram in=$CAT1 x=Total_flux xlog=True ylabel=



## Catalogue Matching using coordinates
topcat -stilts tskymatch2 in1=$CAT1 in2=$CAT2 out=$CAT3 \
    ra1=RA dec1=DEC ra2=right_ascension dec2=declination error=2.0  \
    join=1and2 #join=all1



## matching 3 at once
CAT4=./CATALOGUES/lofar_10sqdeg_edfpos_v4.1_gt5.fits

topcat -stilts tmatchn nin=3 in1=$CAT1 in2=$CAT2 in3=$CAT4 \
    matcher=sky params=2.0 \
    values1="RA DEC" join1=always \
    values2="right_ascension declination" \
    values3="RA DEC" \
    out=./CATALOGUES/triple_matched.fits 

## Example of matching on a column

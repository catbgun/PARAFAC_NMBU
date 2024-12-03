

```{r smoothing and then peaks}
#SMOOTHING
#eemint4smooth <- eem_smooth(eemint4, n = 4, cores = cores)

#WITH SMOOTHING
#bix <- eem_biological_index(eemint4smooth)
#coble_peaks <- eem_coble_peaks(eemint4smooth)
#fi <- eem_fluorescence_index(eemint4smooth)
#hix <- eem_humification_index(eemint4smooth, scale = TRUE)

#indices_peaks_smooth <- bix %>%
#  full_join(coble_peaks, by = "sample") %>%
#  full_join(fi, by = "sample") %>%
#  full_join(hix, by = "sample")

#df_mergeS <- merge(indices_peaks_smooth, met, by = "sample", 
all.x = TRUE) 

#write.csv(df_mergeS,"C:/Users/CBG/OneDrive - NIVA/1 Projects/NMBU_PARAFAC/Output/NMBU_indices and peaks_SMOOTH.csv", row.names = FALSE)
```


HAVE NOT DONE THE FOLLOWING, SKIPPED TO #11)
# 10)Visually find irregularities manually replac e by NA and interpolate
Found several with signal in low bottom left corner. Contamination? Also a few with strange looking signal, but most likely signal close to zero. 
```{r try to remove "protein-like" contamination}
#This cant be done since it interfers with ssample signal elsewhere
#eem_list3x <- eem_setNA(eem_list3, ex =250:300 , em = 250:350, interpolate = 1)
#eemlist3x <- eem_list3 %>% eem_range(ex = c(250,Inf), em = c(0,580))
```

Have not been done: Remove noise, with focus on 2 region
```{r}
eem_list7x <- eem_list7 %>%
  eem_setNA(sample = "DOM11950", ex = 250:275, em=300:350, interpolate = TRUE) %>%
  eem_setNA(sample = "DOM12758", ex = 250:275, em = 300:375, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM13115", ex = 225:275, em = 350:400, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM12695", ex = 250:295, em = 300:375, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM12117", ex = 250:300, em = 300:350, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM12089", ex = 260:285, em = 300:325, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM12864", ex = 250:275, em = 310:325, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM11552", ex = 250:275, em = 300:350, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM11537", ex = 250:275, em = 300:350, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM11389", ex = 260:280, em = 300:340, interpolate = TRUE)%>%
  eem_setNA(sample = "DOM11378", ex = 260:280, em = 340:360, interpolate = TRUE) #usikker

```

#11) First attempt PARAFAC model
minimum and maximum of numbers of components
```{r}
dim_min <- 3
dim_max <- 7

nstart <- 25 # number of similar models from which best is chosen
maxit = 5000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-6 # tolerance in PARAFAC analysis
```

```{r}
pf1 <- eem_parafac(eemint4, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
eempf_compare(pf1, contour = TRUE)
```
# rescale B and C modes to a maximum fluorescence of 1 for each component
```{r}
pf1n <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf1n, contour = TRUE)
```
- har muligens mest tro på modellen med 5 komponenter

From first look it appears to be 4 clean and likely components. There might be 5 or 6 but then it looks like there is something needing cleaning. Or, try with higher PARAFAC resolution? from the guide online we can select 5 or 6. 

#check for correlation between the components
#for a wide range of DOC values, correlation is likely. Then normalise
#The number indicates which model to look at
- there appear to be correlation between the components so selects normalisation
```{r check for correlation}
eempf_cortable(pf1n[[3]], normalisation = FALSE)
eempf_corplot(pf1n[[3]], progress = FALSE, normalisation = FALSE)
```

#Normalise to reduce correlation between components
#wide range of DOC among samples normally require this step
```{r normalise}
pf2 <- eem_parafac(eemint4, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
# rescale B and C modes
pf2nx <- lapply(pf2, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf2nx, contour = TRUE, normalise=FALSE)
```

#Set normalisation TRUE if you want to see the actual data
```{r}
eempf_cortable(pf2nx[[3]], normalisation = FALSE)
eempf_corplot(pf2nx[[3]], progress = FALSE, normalisation = FALSE)
#eempf_compare(pf2nx, contour = TRUE, normalisation=FALSE)
#eempf_plot_comps(pf2, contour = TRUE, type = 1)
```

#Find and exclude outlier leverages
- plot leverage (nice plot)
- plot leverage, not so nice plot but interactive to select what to exclude
- saved in exclude, can be used to start over again with eem_list_ex <- eem_list %>% eem_exclude(exclude) above
- NOt all potential outliers from leverage plot are potential outliers from the other plot. 
```{r calculate leverage}
cpl <- eempf_leverage(pf2nx[[3]])
#cpl <- eempf_leverage(pf3n[[2]])

eempf_leverage_plot(cpl,qlabel=0.1)
#exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

```

should remove sample obenGF317!
  
  ```{r PLot leverage of potential outliers}
eempf_residuals_plot(pf2nx[[3]], eemint4, residuals_only = TRUE, select = c("obenGF317"), spp = 6, cores = cores, contour = TRUE)

# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("obenGF317")
)

# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eemint4, exclude)

```

#remake model - quick
```{r}
pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
# rescale B and C modes
pf3nx <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf3nx, contour = TRUE, normalise=FALSE)

cpl <- eempf_leverage(pf3nx[[3]])

eempf_leverage_plot(cpl,qlabel=0.1)
#exclude <- eempf_leverage_ident(cpl,qlabel=0.1) #sample obenGF164?
eempf_residuals_plot(pf3nx[[3]], eem_list_ex, residuals_only = TRUE, select = c("obenGF164"), spp = 6, cores = cores, contour = TRUE)

```
beholder obenGF164 for nå! eventuelt revisit

#Remake the PARAFAC MODEL
FORTSETT HER!!!! GIR DE MENING ELLER BØR DET RYDDES OPP MER? F EKS SE TILBAKE PÅ SIGNAL I BLANKENE OG PÅ RESIDUALS. KAN VÆRE BEHOV FOR LITT MER OPPRYDNING. SE NÆRMERE PÅ PRØVENE. 
```{r Re-make PARAFAC model}
ctol <- 10^-8     # decrease tolerance in PARAFAC analysis
nstart = 20        # increase number of random starts
maxit = 10000      # increase number of maximum interations

pf4fc <- eem_parafac(eem_list_ex, comps = 5, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, output="all")
pf4fnc <- lapply(pf4fc, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf4fnc, contour = TRUE, normalise=FALSE)

summary(pf4fnc)
#something on the importance of each component
b <- eempf_varimp(pf4fnc[[1]], eem_list7, cores = 2)
b

#eempf_comp_load_plot(pf4fnc[[1]], contour = TRUE)
ggeem(pf4fnc[[1]], contour = TRUE)
```

```{r CHeking and plotting the model}
#leverages, outliers?
cplx <- eempf_leverage(pf4fnc[[1]])
eempf_leverage_plot(cplx,qlabel=0.1)

#residuals, if samples appear suspect from the leverage check
#eempf_residuals_plot(pf4fnc[[1]], eem_list7, residuals_only = TRUE, select = c("DOM11856","DOM12210"), spp = 6, #cores = cores, contour = TRUE)

#plot model
eempf_comp_load_plot(pf4fnc[[1]], contour = TRUE)

#Check the convergence behaviour of the created models:
eempf_convergence(pf4fnc[[1]])

#plot residuals - working?
eempf_residuals_plot(pf4fnc[[1]], eem_list_ex, select = eem_names(eem_list_ex)[10:14], cores = cores, contour = TRUE) #HER ER DET OGSÅ VERDT Å TA EN EKSTRA TITT!!! SER UT SOM EN DEL RESIDUALS OG KANSKJE COMP5 IKKE ER EKTE?

#Plot the resulting components and loadings
eempf_comp_load_plot(pf4fc[[1]], contour = TRUE)
#eempf_residuals_plot(pf4f[[1]], eem_list6, cores = cores, contour = TRUE)

#SPlit-half analysis
sh <- splithalf(eem_list_ex, 5, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, maxit = maxit, ctol = ctol) #HER MÅ DET OGSÅ VURDERES!
splithalf_plot(sh)

#Tucker’s congruency coefficient 
tcc_sh_table <- splithalf_tcc(sh)
tcc_sh_table

#converging models
eempf_convergence(pf4fc[[1]])
```

#To create output file
```{r}
eempf_openfluor(pf4fn[[1]], file = "220924_100Lakes_NEWUVVisdata.txt")
eempf_openfluor(pf3[[3]], file = "221010_100Lakes_NEWUVVisdata.txt")

eempf_openfluor(pf4fnc[[1]], file = "221219_100Lakes_OLDC.txt")
eempf_openfluor(pf4fc[[1]], file = "221211_100Lakes_OLDC2.txt")
```

```{r}
#To create output file
# 9)depending on instrument used, smoothing can be perforemd prior to peak picking, but should not be done prior to PARAFAC
eem4peaks <- eem_smooth(eem_list4, n = 4)
#PEAK PICKING AND INDICES, make sure to use the interpolated version of the data
library(eemR)
bix <- eem_biological_index(eem_list4)
coble_peaks <- eem_coble_peaks(eem_list4)
fi <- eem_fluorescence_index(eem_list4)
hix <- eem_humification_index(eem_list4, scale = TRUE)

#using smoothed data
bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

head(indices_peaks)

write.csv(indices_peaks, "C:/Users/CBG/R/Projects/100 lakes StarDOM/Output/210629_BlanksAsSamples_indexesSmoothed.csv")
```


# 9)depending on instrument used, smoothing can be perforemd prior to peak picking, but should not be done prior to PARAFAC
eem4peaks <- eem_smooth(eem_list5, n = 4)
#PEAK PICKING AND INDICES, make sure to use the interpolated version of the data
library(eemR)
bix <- eem_biological_index(eem_list5)
coble_peaks <- eem_coble_peaks(eem_list5)
fi <- eem_fluorescence_index(eem_list5)
hix <- eem_humification_index(eem_list5, scale = TRUE)

#using smoothed data
bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

head(indices_peaks)

write.csv(indices_peaks, "C:/Users/CBG/R/Output/200819_ALL_noRamannorm_indexes.csv")

#To add TOC values, the sample ID must be coupled to station id
setwd("C:/Users/CBG/R/Projects/100 lakes StarDOM/")
d2 <- read.table("Provenavn.txt", sep = "\t", header = TRUE)
d2$sample = d2$Prøvenr
head(d2)

newtable <- merge(indices_peaks,d2, by  = "sample") 
head(newtable)
newtable <- newtable[, -c(1)]

#Add TOC data and Abs data from file produced for James
setwd("C:/Users/CBG/R/Projects/100 lakes StarDOM/Absorbency spectra")
d1 <- read.table("1000Lakes_Abs (002)_to James.txt", sep = "\t", header = TRUE)
head(d1)

newtable2 <- merge(newtable,d1, by  = "Stasjonskode") 
head(newtable2)

#set columns to keep
new <- newtable2[,c(1:9, 11, 13, 14:19, 22, 24:31)]
head(new)
write.csv(new, "C:/Users/CBG/R/Output/200819_ALL_noRamannorm_indexes.csv")


#Correlation matrix
library("PerformanceAnalytics")
my_data <- new[, c(2,3,4,5,6,7, 8, 9, 20, 22, 23, 24, 25, 26)]
chart.Correlation(my_data, histogram=TRUE, pch=19)

#look at correlation between two variables
library("ggpubr")
ggscatter(new, x = "TOC", y = "Abs254", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")

nrow(new)
write.csv(Dx, "//niva-of5/osl-userdata$/CBG/Documents/R/Projects/100 lakes StarDOM/Output/200805_EEMindices.csv")

#Absorbance index calculation
slope_parms <- abs_parms(absorbance, cuvl = 1)
head(slope_parms)

#connect with the rest
slope_parms$Prøvenr.x = slope_parms$sample

new2 <- merge(new,slope_parms, by  = "Prøvenr.x") 
head(new2)
write.csv(slope_parms, "//niva-of5/osl-userdata$/CBG/Documents/R/Projects/100 lakes StarDOM/Output/200421_abs_indices.csv")



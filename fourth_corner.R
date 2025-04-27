#You need version of R 4.2.2 to work
#Before apporaching, please load datasheet into the working directory
library(readxl)
library(ade4)
library(tidyr)

# Load all sheets into a list
fourth_corner <- lapply(c("sp", "env", "traits"), function(x) read_excel("fourth_corner.xlsx", sheet = x))

# Assign names to the list elements
names(fourth_corner) <- c("sp", "env", "traits")

# Access individual sheets with the following syntax
dim(fourth_corner$sp)     # Dimensions of 'sp' sheet
dim(fourth_corner$env)    # Dimensions of 'env' sheet
dim(fourth_corner$traits) # Dimensions of 'traits' sheet

fourth_corner$env <- as.data.frame(lapply(fourth_corner$env, function(x) {
  if (is.character(x)) as.factor(x) else x
}))

fourth_corner$traits <- as.data.frame(fourth_corner$traits)

fourth_corner$env$Elevation <- as.numeric(as.character(fourth_corner$env$Elevation))
fourth_corner$env <- fourth_corner$env[, !names(fourth_corner$env)]
fourth_corner$traits$Wingspan <- as.numeric(gsub(",", ".", fourth_corner$traits$Wingspan))

fourth_corner$sp[is.na(fourth_corner$sp)] <- 0
fourth_corner$traits[is.na(fourth_corner$traits)] <- 0
afcL.aravo <- dudi.coa(fourth_corner$sp, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(fourth_corner$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(fourth_corner$traits, row.w = afcL.aravo$cw,
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)

nrepet <- 999
four.comb.aravo <- fourthcorner(fourth_corner$env, fourth_corner$sp,
                                fourth_corner$traits, modeltype = 6, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)

plot(four.comb.aravo, alpha = 0.05, stat = "D2")
plot(four.comb.aravo, x.rlq = rlq.aravo, alpha = 0.05,
     stat = "D2", type = "biplot")

Srlq <- fourthcorner2(fourth_corner$env, fourth_corner$sp,
                      fourth_corner$traits,
                      modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

#TIFF extraction of table
tiff('table.tiff', units="in", width=8, height=10, res=600)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()

#PDF is also possible-can be further changed in Adobe Illustrator
pdf('table.pdf', width=8, height=10)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()

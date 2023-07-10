#set working directroy#
setwd("/Users/takeshitaniguchi/R_wd")

#Preparation#
# From "Plantspp_16S_rrarefied_OTU_Table.xlsx" and "Plantspp_ITS_rrarefied_OTU_Table.xlsx", prepare the csv separated file including the sequences collected in summer or winter and name the file of 16S as "rarefied_16S_five_1000_cleaned_dry.csv" and "rarefied_16S_five_1000_cleaned_wet.csv").  

#follwoing is the analysis for 16S in summer#

### Estimation ###
#Analysis of 16S in summer#
# read data #
spp <- read.csv("rarefied_16S_five_1000_cleaned_dry.csv")

#taxnames<-spp[,ncol(spp)]
#names(taxnames)<-rownames(spp)

spp <- spp[,-ncol(spp)]

spp <- subset(spp, rowSums(spp) > 0)

# confirm data structure #
head(spp)
tail(spp)
str(spp)

# get col and row names #
Colnames <- spp$OTUID
Rownames <- colnames(spp[2:ncol(spp)])

# derive only element #
spp <- spp[1:nrow(spp), 2:ncol(spp)]

# convert from dataframe to matrix #
spp <- matrix(as.matrix(spp), nrow(spp), ncol(spp))

# transpose matrix#
spp <- t(spp)

# assign col and row names #
colnames(spp) <- Colnames
rownames(spp) <- Rownames

# library packages #
library("minpack.lm")
library("Hmisc")

# calculate the average number of individuals per community (number of sequences per sumple?) #
N <- mean(apply(spp, 1, sum))

# calculate the relative abundance of each OTU across communities (across sample places) #
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

# calculate the occurence frequency of each OTU across communities (across sample places) #
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq !=0]

# combine the relative abundance and occurence frequency data #
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C <- C[!(apply(C, 1, function(y) any(y==0))),]
p <- C[,2]
freq <- C[,3]

# assign each OTU names to p and freq #
names(p) <- C[,1]
names(freq) <- C[,1]

# calcuation of the limit of detection #
d = 1/N

# estimation of migration parameter using Non-linear least squares #
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))

# get prediction values #
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)

# get 95% confidence interval using Wilson score #
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

# calculate goodness-of-fit (R-squared) #
Rsqr <- 1 - sum((freq - freq.pred)^2)/sum((freq - mean(freq))^2)

# arrange results #
B <- cbind(p, freq, freq.pred, pred.ci[,2:3])

# read key stone list data #
#KS <- read.csv("KeyStoneList.csv")
#KSDF <- data.frame(KeyStone=KS$KeyStone)
#rownames(KSDF) <- KS$OTUID

# merge prediction results and keystone list #
#B <- merge(A, KSDF, by=0)
colnames(B) <- c("p", "freq", "freq.pred", "ci.lower", "ci.upper")

# write prediction data #
write.csv(B, "Result_Prediction_16S_dry.csv")

m.fit 
Rsqr

### Figure ###

library(ggplot2)

# set point colour #
Col <- rep(1,nrow(B))
B <- cbind(B, Col)

B[B$freq > B$ci.upper,]$Col <- 2
B[B$freq < B$ci.lower,]$Col <- 3

B[B$Col==1,]$Col <- "Neutral"
B[B$Col==2,]$Col <- "Above prediction"
B[B$Col==3,]$Col <- "Below prediction"

# set symbol #
#B[B$KeyStone==0,]$KeyStone <- "Others"
#B[B$KeyStone==1,]$KeyStone <- "Key stone"

# create prediction line data #
C <- B[!duplicated(B$p),]

# draw a figure #
p <- ggplot(B, aes(x=log10(p), y=freq, colour=Col))
p <- p + geom_point()
p <- p + scale_shape_manual(values=c(2,16))
p <- p + scale_size_manual(values=c(3,1))
p <- p + scale_colour_manual(values=c("darkgreen", "blue", "black"))
p <- p + geom_line(data=C, aes(x=log10(p), y=freq.pred), col="darkgrey", size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.upper), col="darkgrey", lty=2, size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.lower), col="darkgrey", lty=2, size=1.0)
p <- p + annotate("text", label="R^2*' = 0.77'",x=-3.2, y=0.68, size=5, parse=TRUE)
p <- p + theme_bw()
p <- p + theme(
			   axis.text.x = element_text(size=rel(1.5)),
			   axis.title.x = element_text(size=rel(1.5)),
			   axis.text.y = element_text(size=rel(1.5)),
			   axis.title.y = element_text(size=rel(1.5)),
			   panel.grid.minor = element_blank(),
			   panel.grid.major = element_blank(),
			   legend.position=c(0.51,0.72),
			   legend.justification=c(1,0),
			   legend.text=element_text(size=14),
			   legend.background=element_rect(colour="grey")
			   )
p <- p + xlab(expression(paste(log[10], "(Mean relative abundance)")))
p <- p + ylab("Occurence frequency")
p <- p + guides(scale="none",
				size=FALSE,
				pch=guide_legend(title=NULL),
				colour=guide_legend(title=NULL)
				)
p

#Save the image as 4.7x4.3 incc pdf file#

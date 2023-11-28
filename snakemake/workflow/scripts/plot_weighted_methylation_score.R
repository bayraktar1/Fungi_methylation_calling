library(argparse)
library(karyoploteR)
library(arrow)

parser <- ArgumentParser(description= 'Plot data from feather files with karyoplotr')
parser$add_argument('--genes', '-g', help= 'BED file')
parser$add_argument('--repeats', '-r', help= 'BED file')
parser$add_argument('--feather', '-f', help= 'Binary feather file')
parser$add_argument('--chromosomes', '-c', help= 'Chromosomes and their sizes in tsv')
parser$add_argument('--output', '-o', help = 'Name of output file')
parser$add_argument('--label', '-l', help= 'Name of file')
xargs<- parser$parse_args()

genes_bed <- xargs$genes
repeat_bed <- xargs$repeats
feather_file <- xargs$feather
chromosomes <- xargs$chromosomes
label <- xargs$label
output <- xargs$output


repeats <- toGRanges(repeat_bed)
genes <- toGRanges(genes_bed)
feather_data <- read_feather(feather_file)
feather_ranges <- IRanges(start=feather_data$start, end=feather_data$end, score=feather_data$score)
feather <- GRanges(seqnames=feather_data$Chr, ranges=feather_ranges)
minmeth <- min(feather$score)
maxmeth <- max(feather$score)


df <- read.table(file = chromosomes, sep = '\t', header = FALSE)
pp <- getDefaultPlotParams(plot.type = 4)
pp$leftmargin <- 0.15
atmargin <- 0.20
tracks <- 2

svg(paste(output,".svg"))
chr <- toGRanges(data.frame(chr= df[, 1] , start= seq(from=1,to=nrow(df),by=1), end= df[, 2]))

kp <- plotKaryotype(genome=chr, plot.type=4, plot.params=pp)
kpAddMainTitle(kp, paste("Weighted methylation score", label), cex=1)

at <- autotrack(current.track = 1, total.tracks = tracks, margin = atmargin)
kpPlotRegions(kp, data=reduce(genes), col="#AACCFF", layer.margin = 0.01, border=NA, r0=.1, r1=.2)
kpPlotRegions(kp, data=reduce(repeats), col="#FFDCA8", layer.margin = 0.05, borde=NA, r0=.1, r1=.2)

at <- autotrack(current.track = 2, total.tracks = tracks, margin = atmargin)
kpLines(kp, data=feather, y=feather$score, col="black", ymin = minmeth, ymax = maxmeth, r0=at$r0, r1=at$r1)
kpAxis(kp, r0=at$r0, r1=at$r1, ymin = minmeth, ymax = maxmeth, numticks = 2)
kpAddLabels(kp, labels=label, r0=at$r0, r1=at$r1)
dev.off()


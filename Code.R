#To use your own chain file place the file in C:\Users\00090473\AppData\Local\Temp\RtmpiQwo6v/hg19ToHg38.over.chain


library(glue)
chain_file = "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_file_basename = basename(chain_file)
tmp_dir = tempdir()
local_file = glue("{tmp_dir}/{chain_file_basename}")
download.file(chain_file, destfile = local_file)
system(glue("gzip -d -f '{local_file}'"))


library(rtracklayer)
chain = import.chain(gsub("\\.gz$", "", local_file))

library(circlize)
genome_df = read.chromInfo(species = "hg19")$df
genome_gr = GRanges(seqnames = genome_df[, 1], ranges = IRanges(genome_df[, 2]+1, genome_df[, 3]))


library(EnrichedHeatmap)
genome_bins = makeWindows(genome_gr, w = 5000)

to_bins = liftOver(genome_bins, chain)
to_bins

to_bins = reduce(to_bins, min = 500)
n = sapply(start(to_bins), length)
l = n > 0
to_bins = to_bins[l]
genome_bins = genome_bins[l]

df_from = as.data.frame(genome_bins)[, 1:3]

sq = seqnames(to_bins)
pa = PartitioningByEnd(sq)
df_to = data.frame(
    chr = as.vector(unlist(sq))[start(pa)],
    start = sapply(start(to_bins), function(x) x[1]),
    end = sapply(end(to_bins), function(x) x[1])
)

from_chromosomes = paste0("chr", c(1:22, "X", "Y"))
to_chromosomes = paste0("chr", c(1:22, "X", "Y"))

l = df_from[, 1] %in% from_chromosomes & df_to[, 1] %in% to_chromosomes
df_from = df_from[l, ]
df_to = df_to[l, ]

head(df_from)
head(df_to)
nrow(df_from)
nrow(df_to)
from_chromInfo = read.chromInfo(species = "hg19")$df
from_chromInfo = from_chromInfo[from_chromInfo[, 1] %in% from_chromosomes, ]

to_chromInfo = read.chromInfo(species = "hg38")$df
to_chromInfo = to_chromInfo[to_chromInfo[, 1] %in% to_chromosomes, ]

from_chromInfo[ ,1] = paste0("from_", from_chromInfo[, 1])
to_chromInfo[ ,1] = paste0("to_", to_chromInfo[, 1])
chromInfo = rbind(from_chromInfo, to_chromInfo)

df_from[, 1] = paste0("from_", df_from[, 1])
df_to[, 1] = paste0("to_", df_to[, 1])


chromosome.index = c(paste0("from_", from_chromosomes), paste0("to_", rev(to_chromosomes))) 
chromInfo[, 1] = factor(chromInfo[ ,1], levels = chromosome.index)

library(circlize)
circos.par(gap.after = c(rep(1, length(from_chromosomes) - 1), 5, rep(1, length(to_chromosomes) - 1), 5))
circos.genomicInitialize(chromInfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
        gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

circos.track(ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = "grey", track.height = mm_h(1))

highlight.chromosome(paste0("from_", from_chromosomes), col = "red", track.index = 1)
highlight.chromosome(paste0("to_", to_chromosomes), col = "blue", track.index = 1)

ind = sample(nrow(df_from), 10000)
col = rand_color(length(to_chromosomes))
names(col) = paste0("to_", to_chromosomes)
circos.genomicLink(df_from[ind, ], df_to[ind, ], col = add_transparency(col[df_to[ind, 1]], 0.9))

circos.clear()

text(0.9, -0.9, "hg19")
text(0.9, 0.9, "hg38")


# `to` genome is put on the top of the circle, where the x-axis is reversed
gsize = structure(chromInfo[, 3], names = as.vector(chromInfo[, 1]))
df_to[, 2] = gsize[df_to[, 1]] - df_to[, 2]
df_to[, 3] = gsize[df_to[, 1]] - df_to[, 3]
df_to[, 2:3] = df_to[, 3:2]



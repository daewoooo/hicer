getInteractionSet <- function(h5.file=NULL, chromosomes=NULL, bin.size=10000, min.count=1, bsgenome=NULL, fai=NULL) {
  ## Get chromosome lengths
  if (!is.null(bsgenome)) {
    ## Load BSgenome
    if (class(bsgenome) != 'BSgenome') {
      if (is.character(bsgenome)) {
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- as.object(bsgenome) # replacing string by object
      }
    }
    chrom.lengths <- GenomeInfoDb::seqlengths(bsgenome)
  } else if (!is.null(fai)) {
    ## Get contigs/scaffolds names and sizes from fasta index
    fai.tab <- utils::read.table(fai)
    fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
    chrom.lengths <- fai.tab$V2
    names(chrom.lengths) <- fai.tab$V1
  } else {
    stop("Please submit chromosome lengths in a form of BSgenome object or fasta index (.fai) used for Hi-C alignment !!!")
  } 
  ## TODO or get sequence lengths from the original BAM file
  
  ## Get chromosome to export interaction for
  if (!is.null(chromosomes)) {
    chroms2use <- intersect(names(chrom.lengths), chromosomes)
  } else {
    chroms2use <- names(chrom.lengths)
  }  
  
  ## Define object that specify read pair loading parameters
  res.frags <- GenomicRanges::GRanges()
  GenomeInfoDb::seqlevels(res.frags) <- names(chrom.lengths[chroms2use])
  GenomeInfoDb::seqlengths(res.frags) <- chrom.lengths[chroms2use]
  res.frags.param <- diffHic::pairParam(res.frags, restrict = chroms2use)
  
  if (!is.null(h5.file) & file.exists(h5.file)) {
    message("Counting binned interactions [Time required: medium] ...", appendLF=FALSE); ptm <- proc.time()
    interaction.data <- diffHic::squareCounts(files = h5.file, param = res.frags.param, width = bin.size, filter = as.integer(min.count))
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  } else {
    stop("Input file 'h5.file' not defined or doesn't exists!!!")
  }  
  
  return(interaction.data)
}

#' Plot Hi-C contacts for a user defined genomic region.
#'
#' @param h5.file Hierarchical data format (HDF5) to where Hi-C read pairs are stored.
#' @param region A \code{\link{GRanges-class}} object containing genomic region to plots HIC interactions for.
#' @param min.count Minimum number of Hi-C links between any pair of genomic bins.
#' @param runmed.outlier An odd number to calculate running median value to be used as a new max to remove outliers.
#' @param highlight.gr A \code{\link{GRanges-class}} object containing gnomic region(s) to be highlighted in the final plot.
#' @param position A overall layout of the HIC interactions either 'diagonal' or 'horizontal'.
#' @param color.pal A user defined color palette to be used to color the HiC contact matrix.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @param fai A FASTA index to get lengths of genomic sequences.
#' @param title A user defined character string to be used as a title for the final plot
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' @export
#' 
plotContactMatrixRegional <- function(h5.file=NULL, region=NULL, bin.size=10000, min.count=1, runmed.outlier=0, highlight.gr=NULL, position='diagonal', color.pal=NULL, bsgenome=NULL, fai=NULL, title=NULL) {
  ## Helper function
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  ## Get chromosome lengths
  if (!is.null(bsgenome)) {
    ## Load BSgenome
    if (class(bsgenome) != 'BSgenome') {
      if (is.character(bsgenome)) {
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- as.object(bsgenome) # replacing string by object
      }
    }
    chrom.lengths <- GenomeInfoDb::seqlengths(bsgenome)
  } else if (!is.null(fai)) {
    ## Get contigs/scaffolds names and sizes from fasta index
    fai.tab <- utils::read.table(fai)
    fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
    chrom.lengths <- fai.tab$V2
    names(chrom.lengths) <- fai.tab$V1
  } else {
    stop("Please submit chromosome lengths in a form of BSgenome object or fasta index (.fai) used for Hi-C alignment !!!")
  } 
  
  ## TODO or get sequence lengths from the original BAM file
  
  ## Make sure submitted region do not contain unused seqlevels 
  region.seqname <- unique(as.character(GenomeInfoDb::seqnames(region)))
  region <- GenomeInfoDb::keepSeqlevels(region, value = region.seqname)
  if (is.na(GenomeInfoDb::seqlengths(region))) {
    GenomeInfoDb::seqlengths(region) <- chrom.lengths[region.seqname]
  }
  
  ## Define object that specify read pair loading parameters
  res.frags <- GenomicRanges::GRanges()
  GenomeInfoDb::seqlevels(res.frags) <- names(chrom.lengths)
  GenomeInfoDb::seqlengths(res.frags) <- chrom.lengths
  if (!is.na(GenomeInfoDb::seqlengths(region))) {
    ## Discard regions from outsize of region of interest from read-pair contact calculations
    region.discard <- GenomicRanges::gaps(region)
    region.discard <- region.discard[strand(region.discard) == '*']
    res.frags.param <- diffHic::pairParam(res.frags, discard = region.discard, restrict = region.seqname)
  } else {
    res.frags.param <- diffHic::pairParam(res.frags, restrict = region.seqname)
  }
  
  if (!is.null(h5.file) & file.exists(h5.file)) {
    message("Counting binned interactions [Time required: medium] ...", appendLF=FALSE); ptm <- proc.time()
    interaction.data <- diffHic::squareCounts(files = h5.file, param = res.frags.param, width = bin.size, filter = as.integer(min.count))
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  } else {
    stop("Input file 'h5.file' not defined or doesn't exists!!!")
  }  
  
  ## Get genomic bins for pairs of bins
  bins.first <- InteractionSet::anchors(interaction.data, type="first")
  bins.second <- InteractionSet::anchors(interaction.data, type="second")
  ## Report final HiC interaction object
  n.links <- SummarizedExperiment::assay(interaction.data)
  n.links <- n.links[,1]
  
  ## Reorder interactions
  gi <- InteractionSet::GInteractions(bins.first, bins.second)
  gi <- InteractionSet::swapAnchors(gi)
  bins.first <- InteractionSet::anchors(gi, type="first")
  bins.second <- InteractionSet::anchors(gi, type="second")
  anchor1 <- as.numeric(gi@anchor1)
  anchor2 <- as.numeric(gi@anchor2)
  
  ## Construct data.frame used for plotting
  link.df <- data.frame('M1.seqnames'=as.character(GenomicRanges::seqnames(bins.first)),
                        'M1.start'=as.numeric(start(bins.first)),
                        'M1.end'=as.numeric(end(bins.first)),
                        'M2.seqnames'=as.character(GenomicRanges::seqnames(bins.second)),
                        'M2.start'=as.numeric(start(bins.second)),
                        'M2.end'=as.numeric(end(bins.second)),
                        'value'=n.links,
                        'anchor1'=anchor1,
                        'anchor2'=anchor2)
  
  ## Smooth contact matrix values by setting an outlier values based on running median
  if (runmed.outlier > 0) {
    if ((runmed.outlier %% 2) != 0) {
      outlier <- max(runmed(link.df$value, runmed.outlier)) 
      link.df$value[link.df$value > outlier] <- outlier  
    } else {
      warning("Parameter 'runmed.outlier' has to be an odd number!!!")
    }
  }
  
  if (position == 'diagonal') {
    plt <- ggplot2::ggplot(link.df) +
      geom_rect(aes(xmin=M1.start, xmax=M1.end, ymin=M2.start, ymax=M2.end, fill=value)) +
      #scale_fill_gradient(low = "white", high = "black", trans = 'log', name="# of Interactions (log)") +
      xlab("Chromosome/scaffold name") +
      ylab("Genomic position (Mbp)") +
      theme_bw() +
      theme(aspect.ratio=1) + #Make sure plot is always a square
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    if (!is.null(highlight.gr)) {
      highlight.range.df <- data.frame('x' = start(highlight.gr),
                                       'y' = start(highlight.gr),
                                       'xend' = end(highlight.gr),
                                       'yend' = end(highlight.gr),
                                       'group' = 1:length(highlight.gr))
      highlight.point.df <- data.frame('x' = start(highlight.gr),
                                       'y' = end(highlight.gr),
                                       'group' = 1:length(highlight.gr))
      
      plt <- plt + geom_segment(data=highlight.range.df, aes(x = x, xend = xend,  y =y, yend = yend), color='orange',  size=1) + 
        geom_point(data=highlight.point.df, aes(x = x, y = y), shape=3, color='orange',  size=1, stroke=1)
    }
  } else if (position == 'horizonal') {
    linkHorizontalCoords.df <- diagonal2horizonalHiCcontacts(HICcontacts.df = link.df)
    
    xstart <- min(linkHorizontalCoords.df$x)
    xend <- max(linkHorizontalCoords.df$x)
    ylim2 <- max(linkHorizontalCoords.df$y)
    ## Set plotting theme
    my.theme <- theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
    ## Make the plot
    plt <- ggplot2::ggplot(linkHorizontalCoords.df, aes(x = x, y = y)) +
      geom_polygon(aes(fill = value, group = bin)) +
      #scale_fill_gradient(low = "white", high = "black", trans = 'log', name="Links (log)") + 
      coord_fixed( ratio=1, xlim=c(xstart, xend), ylim=c(0, ylim2) ) +
      scale_y_continuous(labels = comma) +
      scale_x_continuous(labels = comma) +
      xlab("Genomic Position (bp)") +
      ylab("") +
      theme_bw() +
      my.theme
    
    if (!is.null(highlight.gr)) {
      highlight.range.df <- data.frame('x' = start(highlight.gr) - (bin.size / 2),
                                       'y' = 0,
                                       'xend' = end(highlight.gr) - (bin.size / 2),
                                       'yend' = 0,
                                       'group' = 1:length(highlight.gr))
      
      plt <- plt + geom_segment(data=highlight.range.df, aes(x = x, xend = xend,  y =y, yend = yend), color='orange',  size=1)
      #plt + geom_rect(data=highlight.range.df, aes(xmin = x, xmax = xend, ymin = 0, ymax = -bin.size), fill='orange')
    } else {
      stop("HiC plot postion/layout can be either 'diagonal' or 'horizontal' !!!")
    }
  }  
  ## Add color palette
  if (!is.null(color.pal)) {
    if (all(areColors(pal))) {
      plt <- plt + scale_fill_gradientn(colours = color.pal, n.breaks=length(color.pal), trans = 'log', name="Links (log)")
    } else {
      warning("Not all colors defined in 'color.pal' are valid, using default color scheme!!!")
      plt <- plt + scale_fill_gradient(low = "white", high = "black", trans = 'log', name="Links (log)")
    }
  }
  ## Add plot title if defined
  if (!is.null(title)) {
    if (is.character(title)) {
      plt <- plt + ggtitle(title) 
    } else {
      warning("Parameter 'title' has to be a character string in order to be added to the plot!!!")
    }
  }
  ## Return final plot
  return(plt)
}


#' Convert Hi-C contact matrix into coordinates suitable for horizontal visualisation of Hi-C contacts.
#'
#' @param HICcontacts.df A \code{data.frame} object containing binned contact matrices. 
#' Required columns: \itemize{
#' \item{"M1.seqnames"}{mate1 bin chromosome id}
#' \item{"M1.start"}{mate1 bin start}
#' \item{"M1.end"}{mate1 bin end}
#' \item{"M2.seqnames"}{mate2 bin chromosome id}
#' \item{"M2.start"}{mate2 bin start}
#' \item{"M2.end"}{mate2 bin end}
#' \item{"value"}{number of Hi-C read pairs connecting mate1 bin and mate2 bin}
#' \item{"anchor1"}{mate1 bin index ordered by genomic position}
#' \item{"anchor2"}{mate2 bin index ordered by genomic position}}
#' @param yHeightRatio Height ratio between y and x axis.
#' @return A \code{data.frame} object.
#' @author David Porubsky
#' @export
#' 
diagonal2horizonalHiCcontacts <- function(HICcontacts.df=NULL, yHeightRatio=1/3) {
  ## Get binsize of used genomic bins
  bin.size <- (HICcontacts.df$M1.end[1] - HICcontacts.df$M1.start[1]) + 1
  ## Add unique bin id to each 
  HICcontacts.df$binID <- paste0('bin', 1:nrow(HICcontacts.df))
  ## Get region boundaries
  xstart <- min(HICcontacts.df$M1.start)
  xend <- max(HICcontacts.df$M2.end)
  region.width <- xend - xstart
  ## Get differences between matrix anchor indices
  bin.diffs <- abs(HICcontacts.df$anchor1 - HICcontacts.df$anchor2)
  ## Get position of bin centers on x-axis
  centers <- rowMeans(
    data.frame(
      HICcontacts.df$M1.start,
      HICcontacts.df$M2.end ) )
  ## Get bin positions on y-axis
  maxys <- 1 + (bin.diffs - 1) * 0.5
  minys <- maxys - 1
  ## Construct data.frame for plotting
  coords.df <- rbind(
    data.frame(
      x=centers - (bin.size / 2),
      y=rowMeans( data.frame(maxys, minys) ),
      value=HICcontacts.df$value,
      bin=HICcontacts.df$binID ),
    data.frame(
      x=centers,
      y=maxys,
      value=HICcontacts.df$value,
      bin=HICcontacts.df$binID ),
    data.frame(
      x=centers + (bin.size / 2),
      y=rowMeans( data.frame(maxys, minys) ),
      value=HICcontacts.df$value,
      bin=HICcontacts.df$binID ),
    data.frame(
      x=centers,
      y=minys,
      value=HICcontacts.df$value,
      bin=HICcontacts.df$binID ) )
  
  ## Get center of each bin on x-axis
  coords.df$x <- coords.df$x - (bin.size / 2)
  ## Blunt ends of x-axis
  xstart <- xstart - (bin.size / 2)
  xend <- xend - (bin.size / 2)
  ## Scale y-axis to bin.size and set limit based on user defined region of x and y size
  coords.df$y <- pmax(0, coords.df$y)
  coords.df$y <- coords.df$y * bin.size
  ylim2 <- region.width * yHeightRatio
  ## Return final data.frame
  return(coords.df)  
}


#' Extract 'bowtie' shape Hi-C contacts around the inversion breakpoints.
#' 
#' This function takes as an input A \code{\link{InteractionSet-class}} object and extract binned Hi-C interactions
#' around user defined genomic range (inversions range) based on the number of bins to lookup downstream and upstream
#' from the defined genomic range (inversion) breakpoints.
#'
#' @param gi.obj Hierarchical data format (HDF5) to where Hi-C read pairs are stored.
#' @param inversion.gr A \code{\link{GRanges-class}} object containing (likely) inverted genomic region.
#' @param n.bin.lookup A number of bins of certain size to look downstream and upstream from the user defined 'inversion.gr'.
#' @param blacklist.gr A \code{\link{GRanges-class}} object containing gnomic region(s) to be highlighted in the final plot.
#' 
getBowtieContacts <- function(gi.obj=NULL, inversion.gr=NULL, n.bin.lookup=2, blacklist.gr=NULL) {
  ## Get bin.size used in GInteractions object
  bin.size <- unique(width(InteractionSet::anchors(gi.obj, type="first")))
  
  ## Get inversion region boundaries with and without highly identical SDs
  if (class(blacklist.gr) == 'GRanges') {
    inversion.inner.gr <- primatR::subtractRegions(gr = inversion.gr, remove.gr = blacklist.gr, mode = 'flanks')
    mask.flanks <- IRanges::subsetByOverlaps(blacklist.gr, inversion.gr, maxgap = 1000)
    inversion.outer.gr <- GenomicRanges::reduce(c(inversion.gr, mask.flanks))
    ## Extent inversion breakpoint one bin.size on either side of inverted breakpoint
    inversion.outer.gr <- GenomicRanges::resize(inversion.outer.gr, width = width(inversion.outer.gr) + (bin.size * 2), fix = 'center')
  } else {
    inversion.inner.gr <- inversion.gr
    ## Extent inversion breakpoint one bin.size on either side of inverted breakpoint
    inversion.outer.gr <- GenomicRanges::resize(inversion.gr, width = width(inversion.gr) + (bin.size * 2), fix = 'center')
  }  
  
  ## Get breakpoints of both inner and outer inverted region ##
  inversion.inner.breaks.gr <- primatR::getRegionBoundaries(inversion.inner.gr)
  ## Adjust inversion breakpoint to the closest bin boundary
  #start(inversion.inner.breaks.gr) <- floor(start(inversion.inner.breaks.gr) / bin.size) * bin.size
  #end(inversion.inner.breaks.gr) <- start(inversion.inner.breaks.gr)
  bin.pos <- round(start(inversion.inner.breaks.gr) / bin.size) * bin.size
  inversion.inner.breaks.gr <- GenomicRanges::GRanges(seqnames=seqnames(inversion.inner.breaks.gr), ranges = IRanges(start=bin.pos, end=bin.pos))
  
  inversion.outer.breaks.gr <- primatR::getRegionBoundaries(inversion.outer.gr)
  ## Adjust inversion breakpoint to the closest bin boundary
  #start(inversion.outer.breaks.gr) <- floor(start(inversion.outer.breaks.gr) / bin.size) * bin.size
  #end(inversion.outer.breaks.gr) <- start(inversion.outer.breaks.gr)
  bin.pos <- round(start(inversion.outer.breaks.gr) / bin.size) * bin.size
  inversion.outer.breaks.gr <- GenomicRanges::GRanges(seqnames=seqnames(inversion.outer.breaks.gr), ranges = IRanges(start=bin.pos, end=bin.pos))
  
  ## Extract single genomic positions
  inner.break1 <- inversion.inner.breaks.gr[1]
  start(inner.break1) <- start(inner.break1) + 1
  inner.break2 <- inversion.inner.breaks.gr[2]
  outer.break1 <- inversion.outer.breaks.gr[1]
  outer.break2 <- inversion.outer.breaks.gr[2]
  start(outer.break2) <- start(outer.break2) + 1
  
  ## Extend inversions boundary by a certain number of bins
  inner1 <- GenomicRanges::resize(inner.break1, width = (bin.size * n.bin.lookup), fix = 'start')
  inner2 <- GenomicRanges::resize(inner.break2, width = (bin.size * n.bin.lookup), fix = 'end')
  outer1 <- GenomicRanges::resize(outer.break1, width = (bin.size * n.bin.lookup), fix = 'end')
  outer2 <- GenomicRanges::resize(outer.break2, width = (bin.size * n.bin.lookup), fix = 'start')
  
  ## Define bin indices delineating inner and outer region of inversion breakpoint
  right.square <- InteractionSet::linkOverlaps(gi.obj, subject1 = inner1, subject2 = outer2, ignore.strand=TRUE, use.region="both")
  left.square <- InteractionSet::linkOverlaps(gi.obj, subject1 = inner2, subject2 = outer1, ignore.strand=TRUE, use.region="both")
  ## Get interaction counts from the whole interaction object
  n.links <- SummarizedExperiment::assay(gi.obj)
  n.links <- n.links[,1]
  bowtie.links <- sum( c(n.links[right.square$query], n.links[left.square$query]) )
  return(bowtie.links)
}


validateInversionHic <- function(gi.obj=NULL, inversion.gr=NULL, n.bin.lookup=3, blacklist.gr=NULL, n.perm=100, bsgenome=NULL) {
  ## Load BSgenome
  if (!is.null(bsgenome)) {
    if (class(bsgenome) != 'BSgenome') {
      if (is.character(bsgenome)) {
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- as.object(bsgenome) # replacing string by object
      }
    }
  }
  ## Get observed contact around inversion breakpoints
  observed.contacts <- getBowtieContacts(gi.obj = gi.obj, inversion.gr = inversion.gr, n.bin.lookup = n.bin.lookup, blacklist.gr = blacklist.gr)
  
  ## Get observed contact around inversion breakpoints
  rand.counts <- list()
  for (i in 1:n.perm) {
    rand.gr <- primatR::randomizeRanges(gr = inversion.gr, bsgenome = bsgenome)
    rand.count <- getBowtieContacts(gi.obj = gi.obj, inversion.gr = rand.gr, n.bin.lookup = n.bin.lookup, blacklist.gr = blacklist.gr)
    rand.counts[[i]] <- rand.count
  }
  permuted.contacts <- unlist(rand.counts, use.names = FALSE)
  
  ## Test if there is significantly more contact around inversion then at permuted contacts
  wilcox.p.val <- wilcox.test(observed.contacts, permuted.contacts, paired = FALSE, alternative = 'greater')
  
  ## Return observed and permuted values
  results <- list('observed' = observed.contacts,
                  'permuted' = permuted.contacts,
                  'stat.test' = wilcox.p.val)
  return(results)  
}
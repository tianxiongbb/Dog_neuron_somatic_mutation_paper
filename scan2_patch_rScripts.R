#  File share/R/nspackloader.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.r-project.org/Licenses/

local({
    info <- loadingNamespaceInfo()
    pkg <- info$pkgname
    ns <- .getNamespace(as.name(pkg))
    if (is.null(ns))
        stop("cannot find namespace environment for ", pkg, domain = NA);
    dbbase <- file.path(info$libname, pkg, "R", pkg)
    lazyLoad(dbbase, ns, filter = function(n) n != ".__NAMESPACE__.")
})

compute.mutburden <- function (object, gbp.per.genome = get.gbp.by.genome(object), quiet = FALSE)
{
    check.slots(object, c("call.mutations", "depth.profile"))
    muttypes <- c("snv", "indel")
    object@mutburden <- setNames(lapply(muttypes, function(mt) {
        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
        sfp <- object@static.filter.params[[mt]]
        g <- object@gatk[resampled.training.site == TRUE & muttype == mt]
        dptab <- object@depth.profile$dptab
        dptab <- dptab[1:min(max(g$dp) + 1, nrow(dptab)), ]
        if (nrow(g) < 100) {
            warning(paste("only", nrow(g), "resampled germline",
                mt, "sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites"))
            ret <- data.frame(ncalls = NA, callable.sens = NA,
                callable.bp = NA)[c(1, 1, 1), ]
        }
        else {
            s <- object@gatk[pass == TRUE & muttype == mt]
            q = 4
            qbreaks <- quantile(g$dp, prob = 0:q/q)
            # prevent the error that equal values are observed in qbreaks
            for(i in 5:2){if(qbreaks[i]==qbreaks[i-1]){qbreaks[i]=qbreaks[i]+i/10}}
            s$dpq <- cut(s$dp, qbreaks, include.lowest = T, labels = F)
            s$dpq[s$dpq == 3] <- 2
            g$dpq <- cut(g$dp, qbreaks, include.lowest = T, labels = F)
            g$dpq[g$dpq == 3] <- 2
            rowqs <- cut(0:(nrow(dptab) - 1), qbreaks, include.lowest = T, labels = F)
            rowqs[rowqs == 3] <- 2
            qstouse <- c(1, 2, 4)
            s <- s[dpq %in% qstouse]
            g <- g[dpq %in% qstouse]
            ret <- data.frame(ncalls = sapply(qstouse, function(q) sum(s[dpq == q]$pass, na.rm = TRUE)), callable.sens = sapply(qstouse, function(q) mean(g[bulk.dp >= sfp$min.bulk.dp & dpq == q]$resampled.training.pass, na.rm = TRUE)), callable.bp = sapply(split(dptab[, -(1:sfp$min.bulk.dp)], rowqs), sum))
        }
        ret$callable.burden <- ret$ncalls/ret$callable.sens
        ret$rate.per.gb <- ret$callable.burden/ret$callable.bp *
            1e+09/2
        ret$burden <- ret$rate.per.gb * gbp.per.genome
        ret$somatic.sens <- ret$ncalls/ret$burden
        ret$pre.genotyping.burden <- pre.geno.burden
        ret
    }), muttypes)
    object
}

genome.string.to.seqinfo.object <- function (genome = c("hs37d5", "hg38", "mm10", "canFam3"))
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        # return(GenomeInfoDb::Seqinfo(genome = "GRCh37.p13"))
        # NCBI ftp is unconnectable, so use hg19 from UCSC
        t <- GenomeInfoDb::Seqinfo(genome = "hg19")
	seqlevels(t)[1:22] <- 1:22
        return(t)
    }
    else if (genome == "hg38") {
        return(GenomeInfoDb::Seqinfo(genome = "hg38"))
    }
    else if (genome == "mm10") {
        return(GenomeInfoDb::Seqinfo(genome = "mm10"))
    }
    else if (genome == "canFam3") {
        return(GenomeInfoDb::Seqinfo(genome = "canFam3"))
    }
}

genome.string.to.tiling <- function (genome = c("hs37d5", "hg38", "mm10", "canFam3"), tilewidth = 1e+07,
    group = c("auto", "sex", "circular", "all"))
{
    genome <- match.arg(genome)
    group <- match.arg(group)
    if (genome == "hs37d5") {
        species <- "Homo_sapiens"
    }
    else if (genome == "hg38") {
        species <- "Homo_sapiens"
    }
    else if (genome == "mm10") {
        species <- "Mus_musculus"
    }
    else if (genome == "canFam3") {
        species <- "Canis_lupus_familiaris"
    }
    else {
        stop("unsupported genoome string")
    }
    sqi <- genome.string.to.seqinfo.object(genome)
    print("here!!!")
sqi
    # seqinfo is retrieved from UCSC
    chroms.to.tile <- GenomeInfoDb::extractSeqlevelsByGroup(species = species,
	style = seqlevelsStyle(sqi)[1], group = group)
    #chroms.to.tile <- GenomeInfoDb::extractSeqlevelsByGroup(species = species,
    #    style = "NCBI", group = group)
    grs <- GenomicRanges::tileGenome(seqlengths = sqi[chroms.to.tile],
        tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
    grs
}

genome.string.to.bsgenome.object <- function (genome = c("hs37d5", "hg38", "mm10", "canFam3"))
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        require(BSgenome.Hsapiens.1000genomes.hs37d5)
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    }
    else if (genome == "hg38") {
        require(BSgenome.Hsapiens.UCSC.hg38)
        genome <- BSgenome.Hsapiens.UCSC.hg38
    }
    else if (genome == "mm10") {
        require(BSgenome.Mmusculus.UCSC.mm10)
        genome <- BSgenome.Mmusculus.UCSC.mm10
    }
    else if (genome == "canFam3") {
        require(BSgenome.Cfamiliaris.UCSC.canFam3)
        genome <- BSgenome.Cfamiliaris.UCSC.canFam3
    }
    genome
}

genome.to.spmgr.format=c("GRCh37","GRCh38","GRCm37","dog")
names(genome.to.spmgr.format)=c("hs37d5","hg38","mm9","canFam3")

make.scan <- function (single.cell, bulk, genome = c("hs37d5", "hg38", "mm10", "canFam3"),
    region = NULL)
{
    genome <- match.arg(genome)
    new("SCAN2", single.cell = single.cell, bulk = bulk, genome.string = genome,
        genome.seqinfo = genome.string.to.seqinfo.object(genome),
        region = region, gatk = NULL, ab.fits = NULL, ab.estimates = NULL,
        mut.models = NULL, cigar.data = NULL, excess.cigar.scores = NULL,
        static.filter.params = NULL, fdr.prior.data = NULL, fdr = NULL,
        call.mutations = NULL, mutburden = NULL, mutsig.rescue = NULL)
}

# following function is specific to dog
make.integrated.table <- function (mmq60.tab, mmq1.tab, phased.vcf, bulk.sample, sc.samples,
    genome, snv.min.bulk.dp, indel.min.bulk.dp, snv.max.bulk.alt = 0,
    snv.max.bulk.af = 0, indel.max.bulk.alt = 0, indel.max.bulk.af = 0,
    panel = NULL, grs = tileGenome(seqlengths = genome.string.to.seqinfo.object("canFam3")[seqlevels(genome.string.to.seqinfo.object("canFam3"))[1:38]],
    tilewidth = 1e+07, cut.last.tile.in.chrom = TRUE), legacy = FALSE, quiet = TRUE, report.mem = FALSE)
{
    cat("Starting integrated table pipeline on", length(grs),
        "chunks.\n")
    cat("Parallelizing with", future::nbrOfWorkers(), "cores.\n")
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(grs))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i, ]
            pc <- perfcheck(paste("read and annotate raw data",
                i), {
                gatk <- read.tabix.data(path = mmq60.tab, region = gr,
                  quiet = quiet, colClasses = list(character = "chr", integer = "pos"))
                sitewide <- gatk[, 1:7]
                samplespecific <- gatk[, -(1:7)]
                annotate.gatk.counts(gatk.meta = sitewide, gatk = samplespecific,
                  bulk.sample = bulk.sample, sc.samples = sc.samples,
                  legacy = legacy, quiet = quiet)
                annotate.gatk(gatk = sitewide, genome.string = genome,
                  add.mutsig = TRUE)
                annotate.gatk.lowmq(sitewide, path = mmq1.tab,
                  bulk = bulk.sample, region = gr, quiet = quiet)
                annotate.gatk.phasing(sitewide, phasing.path = phased.vcf,
                  region = gr, quiet = quiet)
		print("check point!!!")
		print(panel)
		sitewide
                annotate.gatk.panel(sitewide, panel.path = panel,
                  region = gr, quiet = quiet)
		print("check point!!!")
                annotate.gatk.candidate.loci(sitewide, snv.min.bulk.dp = snv.min.bulk.dp,
                  snv.max.bulk.alt = snv.max.bulk.alt, snv.max.bulk.af = snv.max.bulk.af,
                  indel.min.bulk.dp = indel.min.bulk.dp, indel.max.bulk.alt = indel.max.bulk.alt,
                  indel.max.bulk.af = indel.max.bulk.af, mode = ifelse(legacy,
                    "legacy", "new"))
		print("check point!!!")
            }, report.mem = report.mem)
            p(class = "sticky", amount = 1, pc)
            cbind(sitewide, samplespecific)
        })
    })
    gatk <- rbindlist(xs)
    if (nrow(gatk[somatic.candidate == TRUE]) == 0)
        stop("0 somatic candidates detected. SCAN2 requires somatic candidates to have 0 supporting reads in the matched bulk - perhaps your bulk is too closely related to your single cells?")
    resampling.details <- gatk.resample.phased.sites(gatk)
    list(gatk = gatk, resampling.details = resampling.details)
}

classify.muts <- function (df, genome.string, spectype = "SNV", sample.name = "dummy",
    save.plot = F, auto.delete = T, verbose = FALSE)
{
    if (nrow(df) == 0)
        return(df)
    recognized.spectypes <- c("SNV", "ID")
    if (!(spectype %in% recognized.spectypes))
        stop(sprintf("unrecognized spectype '%s', currently only supporting %s",
            spectype, paste("\"", recognized.spectypes, "\"",
                collapse = ", ")))
    require(SigProfilerMatrixGeneratorR)
    spmgd <- tempfile()
    if (file.exists(spmgd))
        stop(paste("temporary directory", spmgd, "already exists"))
    dir.create(spmgd, recursive = TRUE)
    out.file <- paste0(spmgd, "/", sample.name, ".vcf")
    f <- file(out.file, "w")
    vcf.header <- c("##fileformat=VCFv4.0", "##source=SCAN2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
            "FILTER", "INFO", "FORMAT", sample.name), collapse = "\t"))
    writeLines(vcf.header, con = f)
    s <- df
    mutid <- paste(s$chr, s$pos, s$refnt, s$altnt)
    dupmut <- duplicated(mutid)
    if (verbose)
        cat("Removing", sum(dupmut), "/", nrow(s), "duplicated mutations before annotating\n")
    s <- s[!dupmut, ]
    old.opt <- options("scipen")$scipen
    options(scipen = 10000)
    writeLines(paste(s$chr, s$pos, ".", s$refnt, s$altnt, ".",
        "PASS", ".", "GT", "0/1", sep = "\t"), con = f)
    close(f)
    options(scipen = old.opt)
    if (verbose) {
        mat <- SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(sample.name,
            genome.string, spmgd, seqInfo = TRUE, plot = save.plot)
    }
    else {
        reticulate::py_capture_output(mat <- SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(sample.name,
            genome.string, spmgd, seqInfo = TRUE, plot = save.plot))
    }
    # exclude Y as it is not in the dog reference
    annot.files <- paste0(spmgd, "/output/vcf_files/", spectype,
        "/", c(1:38, "X"), "_seqinfo.txt")
    if (spectype == "ID") {
        colclasses <- c(V2 = "character", V5 = "character", V6 = "character")
    }
    else if (spectype == "SNV") {
        colclasses <- c(V2 = "character")
    }
    annots <- do.call(rbind, lapply(annot.files, function(f) {
        tryCatch(x <- read.table(f, header = F, stringsAsFactors = FALSE,
            colClasses = colclasses), error = function(e) NULL)
    }))
    if (substr(s$chr[1], 1, 3) == "chr")
        annots[[2]] <- paste0("chr", annots[[2]])
    if (spectype == "ID") {
        colnames(annots) <- c("sample", "chr", "pos", "iclass",
            "refnt", "altnt", "unknown")
        newdf <- plyr::join(df, annots[2:6], by = colnames(annots)[-c(1,
            4, 7)])
    }
    else if (spectype == "SNV") {
        colnames(annots) <- c("sample", "chr", "pos", "iclass",
            "unknown")
        newdf <- plyr::join(df, annots[2:4], by = colnames(annots)[-c(1,
            4, 5)])
    }
    if (save.plot) {
        plotfiles <- list.files(paste0(spmgd, "/output/plots/"),
            full.names = T)
        file.copy(plotfiles, ".")
    }
    if (!all(df$chr == newdf$chr))
        stop("df and newdf do not perfectly correspond: df$chr != newdf$chr")
    if (!all(df$pos == newdf$pos))
        stop("df and newdf do not perfectly correspond: df$pos != newdf$pos")
    if (auto.delete)
        unlink(spmgd, recursive = TRUE)
    df$muttype <- newdf$iclass
    df
}

get.gbp.by.genome <- function (object)
{
    if (object@genome.string == "hs37d5") {
        total <- 3137161264
        chrx <- 155270560
        chry <- 59373566
        chrm <- 16571
        return(total - chrx - chry - chrm) * 2/1e+09
    }
    else if (object@genome.string == "hg38") {
        total <- 3209286105
        chrx <- 156040895
        chry <- 57227415
        chrm <- 16569
        return(total - chrx - chry - chrm) * 2/1e+09
    }
    else if (object@genome.string == "mm10") {
        total <- 2730871774
        chrx <- 171031299
        chry <- 91744698
        chrm <- 16299
        return(total - chrx - chry - chrm) * 2/1e+09
    }
    else if (object@genome.string == "canFam3") {
        total <- 2327633984
        chrx <- 123869142
        chry <- 0
        chrm <- 0
        return(total - chrx - chry - chrm) * 2/1e+09
    }
    else {
        warning(paste("gbp not yet implemented for genome", object@genome.string))
        warning("the mutation burden for this analysis is a placeholder!")
        warning("DO NOT USE!")
    }
}

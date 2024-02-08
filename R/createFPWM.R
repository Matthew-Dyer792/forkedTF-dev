#' A function to generate an FPWM class object.
#'
#' This function creates the FPWM object and performs the fork in the motif.
#' @param mainTF [character] character with the name of the main TF.
#' @param partners [character list or vector] List or character vector with the names of the partner TFs.
#' @param cell [character] with the name of the main TF.
#' @param mainTF_MMID [character] with the name of the main TF.
#' @param partners_MMID [character vector] with the name of the main TF.
#' @param forkPosition [numeric] defines the postion in the matrix where the motif will be forked.
#' @param probabilityMatrix [logical] whether the function should return a frequency matrix or probability matrix (Default FALSE).
#' @param scaleFrequencyCounts [logical] whether the count matrix should have equal rowSums across all the rows (Default FALSE).
#' @param flipMatrix [logical] whether to apply reverse complement in case the core motif is after the forkPosition (Default FALSE).
#' @param local_db_path [character] The complete path to the SQLite
#' implementation of TFregulomeR database available at
#' "https://methmotif.org/API_TFregulomeR/downloads/"
#' @examples
#' fpwm <- createFPWM(mainTF ="CEBPB", partners = c("ATF4","ATF7","JUND","CEBPD"), cell = "K562", forkPosition = 5)
#' @return returns a FPWM class object that can be used to plot or write in transfact format.
#' @export
createFPWM <- function(mainTF = NULL,
                       partners = NULL,
                       cell = NULL,
                       mainTF_MMID = NULL,
                       partners_MMID = NULL,
                       forkPosition = NULL,
                       probabilityMatrix = FALSE,
                       scaleFrequencyCounts = FALSE,
                       flipMatrix = FALSE,
                       local_db_path = NULL) {
  tfnames <- FALSE
  tfIDs <- FALSE

  # Check if inputs are compatible
  if (!is.null(mainTF) && !is.null(partners) && !is.null(cell) && !is.null(forkPosition)) {
    tfnames <- TRUE
  }
  if (!is.null(mainTF_MMID) && !is.null(partners_MMID) && !is.null(forkPosition)) {
    tfIDs <- TRUE
  }

  if (tfnames && tfIDs) {
    stop("Incompatible input. Provide 'mainTF', 'partners' and 'cell'. !!OR!! 'mainTF_MMID' and 'partners_MMID'. ")
  }

  if (!(tfnames || tfIDs)) {
    stop("Missing input.\nProvide: 'mainTF', 'partners', 'cell' and 'forkPosition' \nOR\n 'mainTF_MMID', 'partners_MMID' and 'forkPosition'")
  }

  if (probabilityMatrix == TRUE && scaleFrequencyCounts == TRUE) {
    stop("scaleFrequencyCounts only works for count matrices, disable probabilityMatrix.")
  }

  # build api_object
  api_object <- TFregulomeR::.construct_api(
    local_db_path = local_db_path
  )

  if (tfnames) {
    xTF_cell_tissue_name <- suppressMessages(
      TFregulomeR::dataBrowser(
        tf = mainTF,
        cell_tissue_name = cell,
        local_db_path = local_db_path
      )$ID[1]
    )

    partners <- as.list(unlist(partners))
    MMpartners <- partners
    for (i in seq_along(partners)) {
      suppressMessages(
        MMpartners[[i]] <- TFregulomeR::dataBrowser(
          tf = partners[[i]],
          cell_tissue_name = cell,
          local_db_path = local_db_path
        )$ID[1]
      )
    }

    if(is.null(xTF_cell_tissue_name)) {
      stop("Please check the spelling of your mainTF or cell. \nOr there might no be information for this TF-cell combination.")
    }

    if (length(xTF_cell_tissue_name) != 1) {
      stop("More than one record for the combination of TF and cell tissue")
    }

    if(sum(MMpartners %in% partners)) {
      message(paste("The following partners were not found in combination with the cell:", unlist(partners[which(MMpartners %in% partners)])))
      stop("\n")
    }
    peak_id_y_list <- MMpartners
    peak_id_x <- xTF_cell_tissue_name
  }

  # convert probabilityMatrix to motif_type
  if (probabilityMatrix == TRUE) {
    motif_type <- "MEME"
  } else {
    motif_type <- "TRANSFAC"
  }

  if (tfIDs) {
    partners_MMID <- as.list(unlist(partners_MMID))
    peak_id_y_list <- partners_MMID
    peak_id_x <- mainTF_MMID
  }

  if (length(peak_id_y_list) < 2) {
    stop("The list of partner TF should contain 2 or more elements.")
  }

  Motif <- TFregulomeR::intersectPeakMatrix(
    peak_id_x = peak_id_x,
    motif_only_for_id_x = TRUE,
    peak_id_y = peak_id_y_list,
    motif_only_for_id_y = TRUE,
    motif_type = motif_type,
    local_db_path = local_db_path
  )
  # get number of positions in motif
  motif_length <- Motif[1, 1][[1]]@MethMotif_x@MMmotif@width

  if (as.numeric(forkPosition) >= motif_length) {
    stop("The forkPosition is larger than the motif length.")
  }

  if (flipMatrix==TRUE) {
    for( i in 1:length(peak_id_y_list) ){
      Motif[1,i][[1]]@MethMotif_x@MMBetaScore <- Motif[1,i][[1]]@MethMotif_x@MMBetaScore[,motif_length:1] # reverse meth info
      colnames(Motif[1,i][[1]]@MethMotif_x@MMBetaScore) <- rev(colnames(Motif[1,i][[1]]@MethMotif_x@MMBetaScore)) # reverse meth info colnames
      Motif[1,i][[1]]@MethMotif_x@MMBetaScore[,1:(motif_length-1)] <- Motif[1,i][[1]]@MethMotif_x@MMBetaScore[,2:motif_length] # shift the methyation to land on the reverse Cs
      Motif[1,i][[1]]@MethMotif_x@MMBetaScore[,motif_length] <- 0 # pad the last row as we don't know the methylation info
      Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[motif_length:1,] # reverse
      tmp_matrix <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix
      tmp_matrix[,'A'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'T'] # complementary
      tmp_matrix[,'T'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'A'] # complementary
      tmp_matrix[,'C'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'G'] # complementary
      tmp_matrix[,'G'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'C'] # complementary
      Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix <- tmp_matrix

    }
    # adjust fork position
    forkPosition <- which(motif_length:1 %in% forkPosition) - 1
  }

  TFregulomeDataovlaplist <- list()
  TFregulomeDataBetalist <- list()
  TFregulomeDatamatrixlist <- list()
  TFregulomeDatanPeaks <- vector()
  TFregulomeDatanSites <- vector()

  for (i in seq_along(peak_id_y_list)) {
    y <- Motif[1, i]
    TFregulomeDataovlaplist[[i]] <- y[[1]]@overlap_percentage_x
    TFregulomeDataBetalist[[i]] <- y[[1]]@MethMotif_x@MMBetaScore
    TFregulomeDatamatrixlist[[i]] <- y[[1]]@MethMotif_x@MMmotif@motif_matrix
    TFregulomeDatanPeaks[i] <- y[[1]]@MethMotif_x@MMmotif@nPeaks
    TFregulomeDatanSites[i] <- y[[1]]@MethMotif_x@MMmotif@nsites
  }

  FPWM <- new("FPWMClassObj")
  FPWM@xid <- peak_id_x
  FPWM <- updateFPWMClassObj(FPWM, id = peak_id_y_list,
                             nSites = TFregulomeDatanSites,
                             nPeaks = TFregulomeDatanPeaks,
                             matrix = TFregulomeDatamatrixlist,
                             betalevel = TFregulomeDataBetalist,
                             score = TFregulomeDataovlaplist,
                             forkPosition = forkPosition)
  FPWM <- MatrixAdder(
    FPWM,  forkPosition, motif_type, flipMatrix, api_object
  )
  FPWM <- BetaAdder(FPWM, forkPosition)
  FPWM <- ConvertToFTRANSFAC(FPWM, probabilityMatrix, scaleFrequencyCounts)

  # for (i in c(1:length(FPWM@betalevel))) {
  for (i in seq_along(FPWM@betalevel)) {
    X <- FPWM@betalevel[[i]]
    FPWM@betalevel[[i]] <- X[, (forkPosition + 1):ncol(FPWM@betalevel[[i]])]
  }
  FPWM <- ModifyBetaFormat(FPWM)


  # order object based on overlapping score
  FPWM_tmp <- FPWM
  objOrder <- order(unlist(FPWM@score), decreasing = TRUE)
  FPWM_tmp@betalevel <- FPWM@betalevel[objOrder]
  FPWM_tmp@id <- FPWM@id[objOrder]
  FPWM_tmp@matrix <- FPWM@matrix[objOrder]
  FPWM_tmp@score <- FPWM@score[objOrder]

  FPWMPO <- FPWM@forked$PO
  from <- min(FPWMPO[duplicated(FPWMPO)])
  to <- max(FPWMPO)
  ix <- cbind(which(FPWMPO %in% from), which(FPWMPO %in% to))

  for (jx in 1:length(objOrder)) {
    from_row_ix <- ix[objOrder[jx], 1] : ix[objOrder[jx], 2]
    to_row_ix <- ix[jx, 1] : ix[jx, 2]
    FPWM_tmp@forked[to_row_ix, 2:5] <- FPWM@forked[from_row_ix, 2:5]
  }

  # add colnames
  FPWM <- FPWM_tmp
  colnames(FPWM@forked) <- c("PO", "A", "C", "G", "T")
  colnames(FPWM@parentbeta) <- c("PO", "number", "meth")
  for (i in 1:length(FPWM@betalevel)) {
    colnames(FPWM@betalevel[[i]]) <- c("PO", "number", "meth")
  }

  return(FPWM)
}

get_peaks <- function(mm_id, api_object) {
  peak_i <- suppressMessages(
    TFregulomeR:::.loadPeaks(
      id = mm_id,
      includeMotifOnly = TRUE,
      api_object = api_object))

  return(peak_i)
}

create_grange <- function(peak_i) {
  bed_i <- GenomicRanges::GRanges(
    peak_i$chr,
    IRanges::IRanges(
      peak_i$start - 99,
      peak_i$end + 100),
    id = peak_i$id)

  return(bed_i)
}

subset_peaks <- function(bed_y, bed_x) {
  # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
  suppressWarnings(bedx_with_bedy <- subsetByOverlaps(bed_x, bed_y))

  return(bedx_with_bedy)
}


MatrixAdder <- function(fpwmObject, forkPosition, motif_type,
                        flipMatrix, api_object) {
  # get the peaks for the provided TFs
  peak_x <- get_peaks(fpwmObject@xid, api_object)
  peaks_y <- lapply(fpwmObject@id, get_peaks, api_object = api_object)

  # convert peaks to granges
  bed_x <- create_grange(peak_x)
  beds_y <- lapply(peaks_y, create_grange)

  # get peak x which intersects with y
  overlapped_beds <- lapply(beds_y, subset_peaks, bed_x = bed_x)

  # merge the granges
  bedx_with_bedy <- do.call(c, overlapped_beds)

  if (!is(api_object, "API") && !validObject(api_object)) {
    stop("Invalid API object!")
  } else {
    # make the request
    request_content_df <- apiRequest(api_object, id = fpwmObject@xid)
  }

  if (!is.null(request_content_df)) {
    motif_seq_path_x <- request_content_df[1, c("TFBS")]
    motif_seq_x <- read.delim(motif_seq_path_x, sep = "\t", header = FALSE)
  }

  #compute motif matrix
  colnames(motif_seq_x) <- c(
    "chr", "start", "end", "strand", "weight", "pvalue", "qvalue", "sequence"
  )
  motif_seq_x$id <- paste0(fpwmObject@xid,"_motif_sequence_", as.vector(rownames(motif_seq_x)))
  motif_seq_x_grange <- GenomicRanges::GRanges(
    motif_seq_x$chr,
    IRanges::IRanges(
      motif_seq_x$start + 1,
      motif_seq_x$end),
    id = motif_seq_x$id,
    pvalue = motif_seq_x$pvalue,
    sequence = motif_seq_x$sequence
  )

  # extract the highest pvalue sequence from each peak
  overlaps <- GenomicRanges::findOverlaps(
    motif_seq_x_grange,
    bedx_with_bedy,
    select = "first")
  mcols(motif_seq_x_grange)$peak_id <- mcols(bedx_with_bedy)$id[overlaps]
  motif_of_peakx_with_peaky_allInfo <- as.data.frame(motif_seq_x_grange)
  # keep only one sequence per peak (taking the first value)
  motif_of_peakx_with_peaky_allInfo <- motif_of_peakx_with_peaky_allInfo %>%
    dplyr::distinct(peak_id, .keep_all = TRUE)

  # transform to matrix
  motif_of_peakx_with_peaky <- formMatrixFromSeq(
    input_sequence = as.vector(motif_of_peakx_with_peaky_allInfo$sequence),
    motif_format = motif_type)

  if (nrow(motif_of_peakx_with_peaky) > 0 && flipMatrix == TRUE) {
    motif_of_peakx_with_peaky <- motif_of_peakx_with_peaky[nrow(motif_of_peakx_with_peaky):1, ] # reverse
    colnames(motif_of_peakx_with_peaky) <- c("T", "G", "C", "A") # reverse compliment
    motif_of_peakx_with_peaky <- motif_of_peakx_with_peaky[, order(colnames(motif_of_peakx_with_peaky))] # maintain column order
  }

  fpwmObject@parentmatrix <- head(
    motif_of_peakx_with_peaky,
    n = fpwmObject@forkPosition)

  return(fpwmObject)
}

BetaAdder <- function(fpwmObject, forkPosition) {
  S <- fpwmObject@betalevel[[1]][, 1:forkPosition]
  for (i in 2:length(fpwmObject@id)) {
      S <- S + fpwmObject@betalevel[[i]][, 1:forkPosition]
      }
  fpwmObject@parentbeta <- S
  return(fpwmObject)
}

ConvertToFTRANSFAC <- function(fpwmObject, probabilityMatrix, scaleFrequencyCounts) {
  Cnumber <- length(fpwmObject@matrix)
  RowNum <- nrow(fpwmObject@matrix[[1]])
  forkPosition <- fpwmObject@forkPosition
  R <- rep(c((forkPosition + 1):RowNum), times = Cnumber)
  Step <- (RowNum - forkPosition)
  df <- data.frame(colnames(c("PO", "A", "C", "G", "T")))
  df[1:forkPosition, "PO"] <- c(1:forkPosition)
  df[(forkPosition + 1):(length(R) + forkPosition), "PO"] <- R
  df[1:forkPosition, 2:5] <- fpwmObject@parentmatrix
  c <- 1

  for (i in seq(forkPosition + 1, dim.data.frame(x = df)[1], Step)) {
    if (scaleFrequencyCounts != TRUE) {
      df[i:((Step + i) - 1), 2:5] <- fpwmObject@matrix[[c]][(forkPosition + 1):RowNum, ]
    } else {
      scale_target <- sum(df[1, 2:5])
      scale_current <- sum(fpwmObject@matrix[[c]][(forkPosition + 1), ])
      scale_factor <- scale_target / scale_current
      scaled_matrix <- round(fpwmObject@matrix[[c]][(forkPosition + 1):RowNum, ] * scale_factor)
      # check for errors
      df[i:((Step + i) - 1), 2:5] <- t(apply(scaled_matrix, 1, function(x) {
        if (sum(x) != scale_target) {
          x[which.max(x)] <- x[which.max(x)] + (scale_target - sum(x))
        }
        return(x)
      }))
    }
    c <- c + 1
  }

  fpwmObject@forked <- df
  return(fpwmObject)
}

ModifyBetaFormat <- function(fpwmObject) {
  BS1 <- fpwmObject@parentbeta
  BS1 <- cbind(c("beta score<10%", "beta score 10-90%", "beta score>90%"), BS1)
  BS1 <- rbind(c("position", c(1:(ncol(BS1) - 1))), BS1)


  M1 <- matrix(nrow = ((ncol(BS1) - 1) * 3), ncol = 3)


  pos1 <- rep(t(BS1[1, 2:ncol(BS1)]), times = 3)
  M1[, 1] <- pos1

  pos1 <- t(BS1[2, 2:ncol(BS1)])
  M1[1:(ncol(BS1) - 1), 2] <- pos1
  pos1 <- t(BS1[3, 2:ncol(BS1)])
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 2] <- pos1
  pos1 <- t(BS1[4, 2:ncol(BS1)])
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 2] <- pos1

  pos1 <- rep(t(BS1[2, 1]), times = (ncol(BS1) - 1))
  M1[1:(ncol(BS1) - 1), 3] <- pos1
  pos1 <- rep(t(BS1[3, 1]), times = (ncol(BS1) - 1))
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 3] <- pos1
  pos1 <- rep(t(BS1[4, 1]), times = (ncol(BS1) - 1))
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 3] <- pos1

  fpwmObject@parentbeta <- M1

  for (i in seq_along(fpwmObject@betalevel)) {

    BS2 <- as.matrix(fpwmObject@betalevel[[i]])
    BS2 <- cbind(c("beta score<10%", "beta score 10-90%", "beta score>90%"), BS2)
    BS2 <- rbind(c("position", c(fpwmObject@forkPosition + 1:(ncol(BS2) - 1))), BS2)
    M2 <- matrix(nrow = ((ncol(BS2) - 1) * 3), ncol = 3)


    pos2 <- rep(t(BS2[1, 2:ncol(BS2)]), times = 3)
    M2[, 1] <- pos2

    pos2 <- t(BS2[2, 2:ncol(BS2)])
    M2[1:(ncol(BS2) - 1), 2] <- pos2
    pos2 <- t(BS2[3, 2:ncol(BS2)])
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 2] <- pos2
    pos2 <- t(BS2[4, 2:ncol(BS2)])
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 2] <- pos2

    pos2 <- rep(t(BS2[2, 1]), times = (ncol(BS2) - 1))
    M2[1:(ncol(BS2) - 1), 3] <- pos2
    pos2 <- rep(t(BS2[3, 1]), times = (ncol(BS2) - 1))
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 3] <- pos2
    pos2 <- rep(t(BS2[4, 1]), times = (ncol(BS2) - 1))
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 3] <- pos2

    fpwmObject@betalevel[[i]] <- M2}

  return(fpwmObject)
}

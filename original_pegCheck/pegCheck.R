#####################################################
# ~ ZPRI: function to check pegRNA designs ~
# from IDT sequence; by checking exactly sequence ordered there cannot be new mistakes added after
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

# v0
# ! assumes spacer is 20 nt

# v1
# checks vs. target sequence
# give target sequence as 5'-3' genomic strand
# any genomic window that includes all of spacer/PBS/RTT
# ! assumes RHA is at least 9 nucleotides
# check last templated nucleotide

# v2
# returns structured list instead of printing via cat()

library(stringr)

# BiocManager::install("msa")
library(msa)


# small function to remove the r from IDT notation ------------------------

rRemover <- function(peg) {
  peg <- gsub('\\n', '', peg)
  peg <- gsub('r', '', peg)
  return(peg)
}


# DNAtoRNA ----------------------------------------------------------------

DNAtoRNA <- function(peg) {
  peg <- gsub('\\n', '', peg)
  peg <- gsub('T', 'U', peg)
  return(peg)
}

# RNAtoDNA ----------------------------------------------------------------

RNAtoDNA <- function(peg) {
  peg <- gsub('\\n', '', peg)
  peg <- gsub('U', 'T', peg)
  return(peg)
}


# function pegCheck -------------------------------------------------------

pegCheck <- function(peg,
                     target) {

  ### make sure target is all in uppercase
  target <- toupper(target)

  # ! assumes 20-nt spacer

  ### constant scaffold sequences
  stdscaffold  <- 'GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGGACCGAGUCGGUCC'
  std2scaffold <- 'GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC'

  ### alternative scaffold
  altscaffoldAAAA <- 'GUUUUAGAGCUAUGCUGGAAACAGCAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC'
  altscaffoldUAAA <- 'GUUUAAGAGCUAUGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC'
  altscaffoldGAAA <- 'GUUUCAGAGCUAUGCUGGAAACAGCAUAGCAAGUUGAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC'

  ### clean input
  peg <- gsub('\n', '', peg)
  peg <- gsub('r', '', peg)

  ### pegRNA-level sequences
  peg_length <- nchar(peg)
  rna_full   <- peg
  dna_full   <- gsub('U', 'T', rna_full)
  idt_full   <- paste0(sapply(strsplit(rna_full, '')[[1]], function(x) paste0('r', x)), collapse = '')

  ### spacer — first 20 nt
  spacer     <- substr(peg, start = 1, stop = 20)
  spacer_rna <- spacer
  spacer_dna <- gsub('U', 'T', spacer)

  # locate spacer in target to determine PAM strand
  if (str_detect(target, spacer_dna)) {
    pamstrand           <- 'forward'
    spacer_anneal_strand <- 'Reverse'
    spacerwin <- as.vector(str_locate(target, spacer_dna))
    pam       <- substr(target, start = spacerwin[2] + 1, stop = spacerwin[2] + 3)
    if (substr(pam, start = 2, stop = 3) != 'GG')
      stop('PAM is not NGG but ', pam)
    pam_strand_label <- 'Forward'

  } else if (str_detect(reversecomplement(target), spacer_dna)) {
    pamstrand           <- 'reverse'
    spacer_anneal_strand <- 'Forward'
    spacerwin <- as.vector(str_locate(reversecomplement(target), spacer_dna))
    pam       <- substr(reversecomplement(target), start = spacerwin[2] + 1, stop = spacerwin[2] + 3)
    if (substr(pam, start = 2, stop = 3) != 'GG')
      stop('PAM is not NGG but ', pam)
    pam_strand_label <- 'Reverse'

  } else {
    stop('Cannot find spacer sequence in target.')
  }

  ### scaffold
  scaffold_seq_note <- NULL   # nucleotide notation shown in Courier New in app; NULL for most
  if (str_detect(peg, stdscaffold)) {
    scaffold_label <- 'Standard scaffold'
    scaffold <- stdscaffold

  } else if (str_detect(peg, std2scaffold)) {
    scaffold_label <- 'Standard2 scaffold'
    scaffold <- std2scaffold

  } else if (str_detect(peg, altscaffoldAAAA)) {
    scaffold_label <- 'Alternative scaffold (Dang 2015): extra loop, no RNApol pause disruption'
    scaffold <- altscaffoldAAAA

  } else if (str_detect(peg, altscaffoldUAAA)) {
    scaffold_label    <- 'Alternative scaffold (Petri 2022): extra loop, pause disruption'
    scaffold_seq_note <- paste0('guu*A*AAA → guu*U*AAA')
    scaffold <- altscaffoldUAAA

  } else if (str_detect(peg, altscaffoldGAAA)) {
    scaffold_label    <- 'Alternative scaffold (Dang 2015): extra loop, pause disruption'
    scaffold_seq_note <- paste0('guu*A*AAA → guu*G*AAA')
    scaffold <- altscaffoldGAAA

  } else {
    stop('Cannot detect scaffold. Checked: standard, standard2, alternativeAAAA, alternativeUAAA, alternativeGAAA.')
  }
  scaffold_rna    <- scaffold
  scaffold_dna    <- gsub('U', 'T', scaffold)
  scaffold_length <- nchar(scaffold)

  # strip spacer and scaffold to get extension sequence
  extseq <- substr(peg, 21, nchar(peg))
  extseq <- gsub(scaffold, '', extseq)

  ### PBS
  beforecut   <- substr(spacer, start = 1, stop = 17)
  beforecutRC <- reversecomplement(beforecut, mode = 'rna')
  PBSstart    <- substr(beforecutRC, start = 1, stop = 5)

  PBSpos <- str_locate(extseq, PBSstart)[1]
  if (is.na(PBSpos))
    stop('Cannot find start of PBS sequence. Looking for: ', PBSstart)

  PBSend <- substr(extseq, start = PBSpos, stop = nchar(extseq))

  # 3'-end protection
  last4              <- substr(PBSend, start = nchar(PBSend) - 3, stop = nchar(PBSend))
  end_prot_detected  <- (last4 == 'UUUU')
  if (end_prot_detected) {
    PBS <- substr(PBSend, start = 1, stop = nchar(PBSend) - 4)
  } else {
    PBS <- PBSend
  }
  PBS_rna    <- PBS
  PBS_dna    <- gsub('U', 'T', PBS)
  PBS_length <- nchar(PBS)

  # locate PBS in target
  if (pamstrand == 'forward') {
    if (!str_detect(reversecomplement(target), PBS_dna))
      stop('PAM is on Forward strand so expecting PBS on Reverse strand but cannot find it.')
    pbs_anneal_strand <- 'Forward'
    cutpos <- str_locate(target, reversecomplement(PBS_dna))[2] + 1

  } else if (pamstrand == 'reverse') {
    if (!str_detect(target, PBS_dna))
      stop('PAM is on Reverse strand so expecting PBS on Forward strand but cannot find it.')
    pbs_anneal_strand <- 'Reverse'
    cutpos <- str_locate(target, PBS_dna)[1] - 1
  }

  ### RTT
  RTT        <- gsub(PBSend, '', extseq)
  RTT_rna    <- RTT
  RTT_dna    <- gsub('U', 'T', RTT)
  RTT_length <- nchar(RTT)
  lastRT     <- substr(RTT_dna, 1, 1)
  last_nt_ok     <- (lastRT != 'C')
  last_nt_genome <- reversecomplement(lastRT)
  RTins      <- reversecomplement(RTT_dna)

  ### RHA anchor and alignment
  if (pamstrand == 'forward') {
    rha7    <- substr(RTT_dna, start = 1, stop = 7)
    rha7g   <- reversecomplement(rha7)
    ligpos  <- as.vector(str_locate(target, rha7g))[2]
    if (ligpos < cutpos)
      stop('Ligation position (', ligpos, ') is before cut position (', cutpos, ').')
    queryFw <- substr(target, cutpos, ligpos)
    RTinsFw <- RTins

  } else if (pamstrand == 'reverse') {
    rha7g   <- substr(RTT_dna, start = 1, stop = 7)
    ligpos  <- as.vector(str_locate(target, rha7g))[1]
    if (ligpos > cutpos)
      stop('Ligation position (', ligpos, ') is after cut position (', cutpos, ').')
    queryFw <- substr(target, ligpos, cutpos)
    RTinsFw <- reversecomplement(RTins)
  }

  aln      <- msa(c(RTinsFw, queryFw), type = 'dna', method = 'Muscle')
  con      <- msaConsensusSequence(aln)
  consplit <- strsplit(con, '')[[1]]
  edipos   <- which(consplit == '?')

  new      <- as.character(unmasked(aln)[1])
  newsplit <- strsplit(new, '')[[1]]
  old      <- as.character(unmasked(aln)[2])
  oldsplit <- strsplit(old, '')[[1]]

  ### collect edits and RHA
  if (pamstrand == 'forward') {
    edits_df <- data.frame(
      pos = paste0('+', edipos),
      old = oldsplit[edipos],
      new = newsplit[edipos],
      stringsAsFactors = FALSE
    )
    rha_length <- nchar(con) - max(edipos)
    rha_seq    <- substr(con, start = max(edipos) + 1, stop = nchar(con))

  } else if (pamstrand == 'reverse') {
    edits_df <- data.frame(
      pos = paste0('−', (nchar(con) - edipos) + 1),
      old = oldsplit[edipos],
      new = newsplit[edipos],
      stringsAsFactors = FALSE
    )
    lastep     <- abs((nchar(con) - min(edipos)) + 1)
    rha_length <- nchar(con) - lastep
    rha_seq    <- substr(con, start = 1, stop = min(edipos) - 1)
  }

  ### return structured data
  invisible(list(
    peg = list(
      length = peg_length,
      idt    = idt_full,
      rna    = rna_full,
      dna    = dna_full
    ),
    spacer = list(
      length        = nchar(spacer),
      rna           = spacer_rna,
      dna           = spacer_dna,
      anneal_strand = spacer_anneal_strand,
      pam_strand    = pam_strand_label,
      pam           = pam
    ),
    scaffold = list(
      label    = scaffold_label,
      seq_note = scaffold_seq_note,
      length   = scaffold_length,
      rna      = scaffold_rna,
      dna      = scaffold_dna
    ),
    pbs = list(
      length        = PBS_length,
      rna           = PBS_rna,
      dna           = PBS_dna,
      anneal_strand = pbs_anneal_strand
    ),
    rtt = list(
      length         = RTT_length,
      rna            = RTT_rna,
      dna            = RTT_dna,
      last_nt        = lastRT,
      last_nt_ok     = last_nt_ok,
      last_nt_genome = last_nt_genome,
      rt_writes      = RTins,
      rha_length     = rha_length,
      rha_seq        = rha_seq
    ),
    edits = edits_df,
    end_protection = list(
      detected = end_prot_detected
    )
  ))
}


# reversecomplement function ----------------------------------------------
# ChatGPT, then edits

reversecomplement <- function(seq,
                              mode = 'dna') {

  if (mode == 'dna') {
    comps <- c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C')
  } else if (mode == 'rna') {
    comps <- c('A' = 'U', 'U' = 'A', 'C' = 'G', 'G' = 'C')
  } else {
    stop('mode can only be "rna" or "dna"')
  }

  seq  <- toupper(seq)
  cseq <- sapply(strsplit(seq, NULL)[[1]], function(base) comps[base])
  paste(rev(cseq), collapse = "")
}


# complement function -----------------------------------------------------
# ChatGPT, then edits

complement <- function(seq,
                       mode = 'dna') {

  if (mode == 'dna') {
    comps <- c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C')
  } else if (mode == 'rna') {
    comps <- c('A' = 'U', 'U' = 'A', 'C' = 'G', 'G' = 'C')
  } else {
    stop('mode can only be "rna" or "dna"')
  }

  seq  <- toupper(seq)
  cseq <- sapply(strsplit(seq, NULL)[[1]], function(base) comps[base])
  paste(cseq, collapse = "")
}

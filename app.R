library(shiny)
library(bslib)
library(stringr)
library(msa)

source("pegCheck_cc.R")


# UI helpers --------------------------------------------------------------

# One labelled sequence block: small grey label above Courier New sequence
seq_block <- function(label, seq) {
  div(class = "seq-block",
      div(class = "seq-label", label),
      tags$code(class = "seq", seq))
}

# ✅ or ❌ row: emoji stays its natural colour, text is always black (or red on fail)
check_row <- function(text, ok = TRUE) {
  icon       <- if (ok) "\U2705" else "\U274C"
  text_style <- if (ok) "" else "color:#dc3545;"
  div(class = "check-row",
      tags$span(icon, style = "margin-right:6px;"),
      tags$span(text, style = text_style))
}

# Section card with title bar and body
section_card <- function(title, ...) {
  div(class = "section-card",
      div(class = "section-title", title),
      div(class = "section-body", ...))
}


# Build full result UI from pegCheck() return value -----------------------

render_result <- function(res) {

  tagList(

    # ── pegRNA ────────────────────────────────────────────────────────────
    section_card("pegRNA",
      div(class = "stat", paste(res$peg$length, "nt")),
      seq_block("IDT notation", res$peg$idt),
      seq_block("RNA notation", res$peg$rna),
      seq_block("DNA notation", res$peg$dna)
    ),

    # ── Spacer ────────────────────────────────────────────────────────────
    section_card("Spacer",
      div(class = "stat", paste(res$spacer$length, "nt")),
      seq_block("RNA notation", res$spacer$rna),
      seq_block("DNA notation", res$spacer$dna),
      div(class = "checks",
          check_row(paste("Protospacer anneals on", res$spacer$anneal_strand, "strand")),
          check_row(paste("PAM is on", res$spacer$pam_strand, "strand and is", res$spacer$pam))
      )
    ),

    # ── Scaffold ──────────────────────────────────────────────────────────
    section_card("Scaffold",
      div(class = "check-row",
          tags$span("\U2705", style = "margin-right:6px;"),
          tags$span(
            res$scaffold$label,
            if (!is.null(res$scaffold$seq_note))
              tagList(" ", tags$code(class = "seq-inline", res$scaffold$seq_note)),
            paste0(" (", res$scaffold$length, " nt)")
          )
      ),
      seq_block("RNA notation", res$scaffold$rna),
      seq_block("DNA notation", res$scaffold$dna)
    ),

    # ── PBS ───────────────────────────────────────────────────────────────
    section_card("PBS",
      div(class = "stat", paste(res$pbs$length, "nt")),
      seq_block("RNA notation", res$pbs$rna),
      seq_block("DNA notation", res$pbs$dna),
      div(class = "checks",
          check_row(paste("PBS anneals on", res$pbs$anneal_strand, "strand"))
      )
    ),

    # ── RTT ───────────────────────────────────────────────────────────────
    section_card("RTT",
      div(class = "stat", paste(res$rtt$length, "nt")),
      seq_block("RNA notation", res$rtt$rna),
      seq_block("DNA notation", res$rtt$dna),
      div(class = "subsection",
          seq_block(paste0("RHA — ", res$rtt$rha_length, " nt"), res$rtt$rha_seq)
      ),
      div(class = "subsection",
          seq_block("RT will write in genome", res$rtt$rt_writes)
      ),
      div(class = "checks",
          check_row(
            paste0("Last nucleotide read by RT is ", res$rtt$last_nt,
                   " — incorporated as ", res$rtt$last_nt_genome, " in the genome",
                   if (!res$rtt$last_nt_ok)
                     " (C not recommended, see Anzalone 2019 Extended Data Fig. 5)"
                   else ""),
            ok = res$rtt$last_nt_ok
          )
      )
    ),

    # ── Edited positions ──────────────────────────────────────────────────
    section_card("Edited positions",
      div(class = "edit-note",
          "Forward genomic strand. Left of nick = −1, right of nick = +1."),
      do.call(tagList, lapply(seq_len(nrow(res$edits)), function(i) {
        row <- res$edits[i, ]
        div(class = "edit-row",
            tags$span(class = "edit-pos", paste0("position ", row$pos)),
            tags$span(":", style = "margin:0 6px; color:#aaa;"),
            tags$code(class = "nuc", row$old),
            tags$span("→", style = "margin:0 8px; color:#888;"),
            tags$code(class = "nuc", row$new))
      }))
    ),

    # ── 3'-end protection ─────────────────────────────────────────────────
    section_card("3′-end protection",
      if (res$end_protection$detected) {
        tagList(
          check_row("3′-end protection detected"),
          seq_block("RNA notation", "UUUU"),
          seq_block("DNA notation", "TTTT")
        )
      } else {
        div(class = "check-row",
            tags$span("No 3′-end protection detected"))
      }
    )

  )
}


# UI ----------------------------------------------------------------------

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    body { background-color: #f5f6f8; }
    .main-card { max-width: 820px; margin: 36px auto 60px; }
    textarea.form-control { font-family: 'Courier New', Courier, monospace; font-size: 0.88em; }

    /* ── section cards ── */
    .section-card {
      border: 1px solid #dee2e6;
      border-radius: 8px;
      margin-bottom: 14px;
      overflow: hidden;
    }
    .section-title {
      background: #f0f4f8;
      padding: 7px 16px;
      font-weight: 700;
      font-size: 0.95em;
      border-bottom: 1px solid #dee2e6;
      letter-spacing: 0.02em;
    }
    .section-body {
      padding: 14px 16px;
    }
    .stat {
      color: #666;
      font-size: 0.88em;
      margin-bottom: 10px;
    }

    /* ── sequence blocks ── */
    .seq-block { margin-bottom: 10px; }
    .seq-label {
      font-size: 0.70em;
      color: #bbb;
      text-transform: uppercase;
      letter-spacing: 0.07em;
      margin-bottom: 1px;
    }
    code.seq {
      font-family: 'Courier New', Courier, monospace !important;
      font-size: 0.9em;
      color: #333;
      word-break: break-all;
      white-space: pre-wrap;
      display: block;
      background: none;
      padding: 0;
      border: none;
    }

    /* ── subsections within RTT ── */
    .subsection {
      border-top: 1px solid #f0f0f0;
      margin-top: 10px;
      padding-top: 10px;
    }

    /* ── check rows ── */
    .checks { margin-top: 10px; }
    .check-row {
      margin: 4px 0;
      font-size: 0.92em;
      display: flex;
      align-items: baseline;
    }

    /* ── edit rows ── */
    .edit-note {
      font-size: 0.80em;
      color: #999;
      margin-bottom: 10px;
    }
    .edit-row {
      display: flex;
      align-items: center;
      margin: 4px 0;
      font-size: 0.92em;
    }
    .edit-pos {
      font-family: inherit;
      font-size: 0.92em;
      color: #333;
      min-width: 110px;
    }
    code.seq-inline {
      font-family: 'Courier New', Courier, monospace;
      font-size: 0.88em;
      color: #333;
      background: none;
      padding: 0;
      border: none;
    }
    code.nuc {
      font-family: 'Courier New', Courier, monospace;
      font-size: 0.92em;
      background: #f4f4f4;
      padding: 1px 6px;
      border-radius: 3px;
      color: #333;
    }
  "))),

  div(class = "main-card",

    card(
      card_header(h4("pegCheck", class = "mb-0 fw-bold")),
      card_body(
        padding = "1.4rem",

        textAreaInput(
          "peg",
          label       = "pegRNA sequence",
          value       = "",
          rows        = 4,
          width       = "100%",
          placeholder = "Paste IDT sequence (rArUrGr… or plain RNA/DNA)…"
        ),

        tags$div(
          class = "mb-3",
          tags$label(HTML("5′–3′ genomic target"), class = "form-label"),
          tags$div(
            "Sequence window of any size, as long as it encompasses spacer, RTT, and PBS.",
            style = "font-size:0.82em; color:#999; margin-bottom:5px;"
          ),
          textAreaInput(
            "target",
            label       = NULL,
            value       = "",
            rows        = 4,
            width       = "100%",
            placeholder = "Paste 5′-3′ DNA sequence…"
          )
        ),

        actionButton(
          "check", "CHECK",
          class = "btn btn-primary w-100",
          style = "font-size:1.1em; letter-spacing:0.06em; padding:10px;"
        )
      )
    ),

    uiOutput("result_panel")
  )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  observeEvent(input$check, {

    peg_val    <- trimws(input$peg)
    target_val <- trimws(input$target)

    if (nchar(peg_val) == 0 || nchar(target_val) == 0) {
      output$result_panel <- renderUI(
        div(class = "alert alert-warning mt-3",
            "Please fill in both fields before clicking CHECK.")
      )
      return()
    }

    result <- tryCatch({
      res <- pegCheck(peg_val, target_val)
      list(ok = TRUE, data = res)
    }, error = function(e) {
      msg <- gsub("^[\t > ]+", "", conditionMessage(e))
      list(ok = FALSE, msg = msg)
    })

    output$result_panel <- renderUI({
      if (result$ok) {
        div(style = "margin-top:16px;", render_result(result$data))
      } else {
        div(
          class = "alert alert-danger mt-3",
          tags$strong("Error — "),
          result$msg
        )
      }
    })
  })
}


shinyApp(ui, server)

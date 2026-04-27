# render_theorem.ps1
#
# Render THEOREM.md to theorem.pdf using pandoc.
#
# Prerequisites:
#   1. pandoc        https://pandoc.org/installing.html
#   2. A LaTeX engine. The default below uses xelatex, available via:
#        - MiKTeX     https://miktex.org/download
#        - TeX Live   https://www.tug.org/texlive/
#
# Fallback if no LaTeX engine is available:
#   Replace the pdf-engine line below with `--pdf-engine=wkhtmltopdf` after
#   installing wkhtmltopdf from https://wkhtmltopdf.org/downloads.html.
#   Math rendering is less polished but the document will produce.
#
# Usage:
#   pwsh ./render_theorem.ps1
# or simply:
#   ./render_theorem.ps1

$ErrorActionPreference = "Stop"

$source = "THEOREM.md"
$output = "theorem.pdf"

if (-not (Get-Command pandoc -ErrorAction SilentlyContinue)) {
    Write-Error "pandoc not found on PATH. Install from https://pandoc.org/installing.html"
    exit 1
}

Write-Host "Rendering $source -> $output ..."
pandoc $source `
    -o $output `
    --pdf-engine=xelatex `
    -V geometry:margin=1in `
    -V mainfont="Times New Roman" `
    -V fontsize=11pt `
    --toc

if ($LASTEXITCODE -eq 0) {
    Write-Host "Wrote $output"
} else {
    Write-Error "pandoc failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

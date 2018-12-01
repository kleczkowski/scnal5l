# -*- cperl -*-
@default_files = ('main.tex');

$pdflatex = 'lualatex --shell-escape --interaction=nonstopmode --halt-on-error %O %S';
$pdf_mode = 1;
$postscript_mode = $dvi_mode = 0;

$pdf_previewer = 'xpdf -remote %R %O %S';
$pdf_update_method = 4;
$pdf_update_command = 'xpdf -remote %R -reload';

;; Ce qui suit est un "manifeste" équivalent à la ligne de commande que vous avez donnée.
;; Vous pouvez le stocker dans un fichier que vous pourrez ensuite passer à n'importe quelle
;; commande 'guix' qui accepte une option '--manifest' (ou '-m').

(specifications->manifest
  (list "r-minimal" "r-renv" "r-rcurl" "r-cairo" "r-stringi" "r-svglite" "r-deseq2" "r-ggplot2" "r-ggrepel" "r-ihw" "r-complexheatmap" "r-tximport" "r-dplyr" "r-tidyr" "r-circlize" "r-venndiagram" "coreutils" "gcc-toolchain"))

# TP Prédiction de gènes


## Résultat

Diagramme de Venn : \
![](/output/jVenn_chart.png)
  
Les gènes correspondent à une sous-séquence des transcripts traduites en protéines par le ribosome. Ils comportent un cadre de lecture consistant en des triplets consécutifs depuis un codon d’initiation ( 'AUG', 'UUG', 'CUG', 'AUU' ou 'GUG') et un codon stop (UAA', 'UAG', ou 'UGA'). Ces codons sont dans le même cadre de lecture ! 
On retrouve en amont du codon d’initiation un motif permettant l'initiation de la traduction via la fixation de la sous-unité 16S de l’ARN ribosomique : AGGAGGUAA appelée séquence de Shine-Dalagarno  [Shine et Dalgarno 1974]. Ce motif n'est pas nécessairement dans le même cadre de lecture que le codon d'initiation et peut être incomplet.

## Sujet

Peu d’organismes bénéficient à ce jour d’une annotation vérifiée expérimentalement. La prédiction des gènes reste donc une tâche importante pour l’annotation automatique des génomes. De multiples logiciels et approches techniques existent pour cette tâche :
https://en.wikipedia.org/wiki/List_of_gene_prediction_software

Nous développerons de notre côté une approche simple pour prédire les gènes des procaryotes basée sur la détection des cadres de lecture et du motif Shine-Dalgarno. L’objectif de ce TP sera de prédire les gènes du génome de référence de Listeria monocytogenes EGD-e (assemblée et séquencée par l’Institut Pasteur), qui présente 2867 gènes:
https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/159/159660%7CListeria%20monocytogenes%20EGD-e/

git clone https://github.com/aghozlane/geneprediction-tpt
Et ajoutez ce travail dans votre dépôt github (à créer en ligne):
git remote add mydepot https://github.com/mylogin/geneprediction-tp
git push -u mydepot master

Dans le dossier geneprediction-tp/data/, vous trouverez:
listeria.fna : génome de Listeria monocytogenes
proteins_665_300274.csv : fiche NCBI des protéines de listeria
prodigal.csv: position des gènes prédites par prodigal
position.csv: position des gènes vérifiée expérimentalement



Vous éditerez le programme Python3 nommé gpred.py dans le dossier gpred/.  Il prendra en argument:
- i fichier fasta avec des séquences nucléotidiques (.fna)
-g la longueur minimale des gènes
-s la distance maximale du codon initiateur où chercher la séquence de Shine-Dalgarno
-d le gap minimal entre deux gènes (Séquence de Shine-Dalgarno indépendant)
-p fichier de sortie des positions brutes des gènes
-o fichier fasta donnant les séquences fasta des gènes prédits

Votre programme utilisera impérativement ces arguments.
Nous utiliserons les espaces en unique marqueurs de l’indentation de votre code. 

Vous utiliserez la librairie pytest pour valider vos développements
pip3 install --user pytest pytest-cov 


La création d’un environnement conda peut résoudre les soucis d’installation (si besoin).
conda create -n gpred python
conda activate gpred 


Vous testerez vos fonctions à l’aide de la commande pytest --cov=gpred -v à exécuter dans le dossier geneprediction-tp/.
En raison de cette contrainte, les noms des fonctions ne seront pas libres. Il sera donc impératif de respecter le nom des fonctions “imposées”, de même que leurs caractéristiques et paramètres. 

Lecture du fichier de séquences /1

read_fasta /1 prend un seul argument correspondant au fichier fasta (fasta_file) et retourne une séquence sous la forme d’un string (sans retour chariot). Par souci de simplicité, le fichier fasta ne contiendra qu’une seule molécule d’ADN. read_fasta s’assurera que les caractères de cette string sont bien en majuscule.

Recherche des codons d’initiation, stop et de la séquence de Shine-Dalgarno /13

La recherche des motifs sera effectuée à l’aide de regex compilées générés par la fonction re.compile. Nous remplacerons l‘uracile par la thymine car nous effectuons cette recherche directement sur le génome ! 
regexp_object = re.compile(“myregex”)
Deux fonctions s’appliquant aux regex_object vont nous intéresser:
regexp_object.search(sequence, [start, [stop]]) : Cherche la première occurrence d’une regex entre une position start et stop de la séquence et retourne un objet: match_object
regexp_object.finditer(sequence, [start,[stop]]): retourne un itérateur de match_object non-chevauchant

Les matchs objects permettent d’identifier la position du match à l’aide des fonctions:
.start(0) ou .end(0)

find_start /4 prend 4 arguments
start_regex: regex object permettant d’identifier un codon start
sequence: Séquence du génome
start : position de début de la recherche
stop : position de fin de la recherche
find_start retourne la position de la première occurrence d’un codon d'initiation dans la zone de recherche et None si rien n’est identifié. S

find_stop /4
stop_regex: regex object permettant d’identifier un codon stop
sequence: Séquence du génome
start : position de début de la recherche
find_stop retourne le premier codon stop se trouvant dans le même cadre de lecture que le codon d'initiation et None si rien n’est identifié

has_shine_dalgarno /5
shine_regex: regex object permettant d’identifier une séquence de Shine-Dalagarno
sequence: Séquence du génome
start: position de début de la recherche
max_shine_dalgarno_distance : Position relative de la séquence de Shine-Dalgarno par rapport au codon d’initiation (valeur obtenue en argument du programme)

has_shine_dalgarno  retourne True si un motif de Shine-Dalgarno a été identifié entre max_shine_dalgarno_distance et -6 nucléotides en amont du codon d’initiation et False sinon. La séquence de Shine-Dalgano ne se trouve pas nécessairement dans le même cadre de lecture que le codon d’initiation.

Algorithme de recherche /4

 predict_genes /4
sequence : Séquence du génome
start_regex:  regex object permettant d’identifier un codon start
stop_regex:  regex object permettant d’identifier un codon stop
shine_regex: regex object permettant d’identifier une séquence de Shine-Dalagarno
min_gene_len: longueur minimale d’un gène (valeur obtenue en argument du programme)
max_shine_dalgarno_distance: Position relative de la séquence de Shine-Dalgarno par rapport au codon d’initiation (valeur obtenue en argument du programme)
min_gap: Distance relative minimale entre 2 gènes (valeur obtenue en argument du programme)

A partir des entrées et fonctions réalisées précédemment, une implémentation de la recherche de gènes peut être réalisée par programmation dynamique: 
position_courante = 0
Tant que longueur_sequence - position_courante >= min_gap
	position_courante = Trouver un codon d’initiation à partir de la position courante
	Si position courante non nulle:
stop = Trouver un codon stop à partir de la nouvelle position courante
Si stop est non nulle:
	Déterminer si le gène répond au critère de longueur minimale
	Si oui:
Déterminer s’il dispose d’une séquence de Shine-Dalgarno en amont de son codon d’initiation
Si oui:
gène probable identifié
position_courante = position dernière lettre codon stop + min_gap
				Sinon:
					position_courante = position_courante + 1
		Sinon:
			position_courante = position_courante + 1
		Sinon:
			position_courante = position_courante + 1
Retourne la liste des positions des gènes prédits

predict_genes retourne la liste des gènes prédits sous la forme d’une liste de listes:
[ [15,80], [250, 440], …]

Attention: La position d’un nucléotide sur une séquence génomique démarre à 1 et termine par la position de la dernière lettre du codon stop !

Programme principal /2

Votre programme principal fera appel aux fonctions précédentes nécessaires pour chercher les gènes dans le sens 5’ puis 3’, puis fera appel à reverse complement pour faire cette recherche dans le sens 3’ vers 5’. Aucune modification de vos fonctions ne doit être réalisée pour effectuer cette seconde partie.
La position des gènes prédits dans le sens 3’ vers 5’ devra être corrigée pour qu’elle apparaisse dans le sens 5’ 3’  ! Les positions ont une valeur croissante uniquement.
Vérifiez la qualité de vos prédictions sur jvenn (en copier coller les prodigal/position/predict_genes):
jvenn (inra.fr)


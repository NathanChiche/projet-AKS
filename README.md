# projet-AKS

Comment déterminer de façon déterministe et en temps polynomial en la taille d'un nombre s'il est premier? L'algorithme AKS répond à cette question et 
nous avons essayé d'implémenter de manière assez efficace cet algorithme en C.

Pour compiler le fichier: gcc akstest-7.c -lmpfr -lgmp -o akstest
Pour l'exécuter: ./akstest

La sortie donnera quelques informations sur le nombre d'entrée et la partie de l'algortihme qui permet de conclure.

Il faut avoir téléchargé la bibliothèque gmp au préalable (le code a été fait pour pouvoir faire usage de grands nombres).

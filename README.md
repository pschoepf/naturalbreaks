# Naturalbreaks
Ports of fisher/jenks natural breaks algorithm with O(k×n×log(n)) complexity, engineered and described by Maarten Hilferink:
http://wiki.objectvision.nl/index.php/Fisher%27s_Natural_Breaks_Classification  

The code provided is a (nearly 1:1, as there is not much left for optimization) port of the corresponding C code that was
originally developed by Maarten Hilferink, © Object Vision BV and which can be found here:
http://wiki.objectvision.nl/index.php/CalcNaturalBreaksCode

My contribution is:
- a basic port to Java (src/main/java)
- a basic port to javascript (src/main/javascript) which allows direct execution in the browser to compute natural looking color scales.

Please see the contained unit tests for usage.

License is provided on the same GNU GPL v3.0 license as the original contribution from Marteen.

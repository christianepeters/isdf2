# Information-Set Decoding (ISD)

Information-Set Decoding is the best known attack against the McEliece Cryptosystem.
The Classic McEliece system (see also 
[https://classic.mceliece.org/](https://classic.mceliece.org/))
uses binary codes, thus codes over the finite field
F<sub>2</sub>.


## Related Publications

* [BLP2008]
Daniel J. Bernstein, Tanja Lange, Christiane Peters.
**Attacking and defending the McEliece cryptosystem.**
In Post-Quantum Cryptography, Lecture Notes in Computer Science,
Vol. 5299, pp. 31–46. Springer-Verlag Berlin Heidelberg, 2008.
[pdf](http://eprint.iacr.org/2008/318)
[doi](http://www.springerlink.com/content/68v69185x478p53g/)
[bibtex](http://dblp.uni-trier.de/rec/bibtex/conf/pqcrypto/BernsteinLP08)
[press release](http://www.hyperelliptic.org/tanja/press/mceliece.html)

* [BLP2011] Daniel J. Bernstein, Tanja Lange, Christiane Peters. **Smaller decoding exponents: ball-collision decoding.**
In CRYPTO 2011, Lecture Notes in Computer Science, Vol. 6841, pp. 743–760. Springer-Verlag Berlin Heidelberg, 2011.
[pdf](http://eprint.iacr.org/2010/585.pdf)
[doi](http://www.springerlink.com/content/9138k05502234348/)
[bibtex](http://dblp.uni-trier.de/rec/bibtex/conf/crypto/BernsteinLP11)


## Iteration count scripts

We use scripts to estimate the complexity of
**information-set decoding attacks for binary codes** and to
determine parameters for the McEliece Cryptosystem.

**Complexity estimations in [BLP2008]**:

* [type1.c](type1.c) and [type3.c](type3.c) are written in C
  using the [MPFI library](https://gforge.inria.fr/projects/mpfi/).

* [isdf2.gp](isdf2.gp) is to be used with the
  [PARI/GP](https://pari.math.u-bordeaux.fr/) computer algebra system.
  The counts are a little less precise.

**Bounds as in [BLP2011]**:

* [blp2011-bound.sage](blp2011-bound.sage) is to be used with the
  open source mathematics software
  [Sage](https://www.sagemath.org/) which is largely Python
  based. 

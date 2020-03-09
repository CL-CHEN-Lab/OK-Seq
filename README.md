# OK-Seq
R package for the analysis of OK-Seq data to study DNA replication fork directionality: from count matrices, RFD calculation to inititation/termination zone calling.

Package: hmmPolarity

Title: Generate the replication fork direcrionality genome-wide and initiation zones calling for OK-Seq

Version: 0.1.0

Author: Yves d'AUBENTON-CARAFA <yves@daubenton.net>, Chun-long CHEN <chunlong.chen@curie.fr>

Maintainer: Yaqun LIU <yaqun.liu@curie.fr>

Description: This R package is served for analysing OK-Seq data. Firstly it could transform data into RFD (replication fork directionality) profiles for a primary visualisation and then by using the HMM package of R (http://www.r-project.org/), it can identify accurately most of the replication initiation zones, termination zones and also the intermediate states.

Depends: R (>= 3.1.0)

Imports:
HMM

License: GNU General Public License v3.0

Encoding: UTF-8

LazyData: true

RoxygenNote: 7.0.2

# Vimtilo
### Vale e IFES, Mapeamento em Trilhos com Intervenções LOcalizadas.


Esse repositório contém o algoritmo desenvolvido na dissertação **ESTIMATIVA DOS ENTALHES DE DORMENTES EM PONTES SEM LASTRO UTILIZANDO SENSOR DE PERFIL A LASER**
dividido em duas pastas.

- O primeira, **Ponte-laboratório-3D**, possui o algoritmo e o banco de dados utilizados para o seu processamento o qual foi adquirido no laboratório do instituto federal do espirito santo campus serra com auxílio de uma impressão em escala 1:10 de uma ponte sem lastro.

- O segundo, **Ponte-Real**, possui o algoritmo e o banco de dados utilizados para o seu processamento o qual foi adquirido em uma estrutura de ponte real durante o período de intervenção da via para manutenção gerenciado pela Vale S.A.

Em ambas pastas o arquivo **main** contém o código principal os outros .m existentes são funções que auxiliam no processamento de dados.

Para realizar corretamente os testes em seu computador é necessário realizar o download dos arquivos e abri-los pelo software matlab executando o arquivo main

Obs: o *Current Folder* do *MatLab* deve estar direcionado para o destino do arquivo main que deseja ser executado, caso contrário o *MatLab* não localizara o banco de informação dos dados da nuvem de pontos


## Software de Aquisição de dados

links para download dos software de aquisição de dados do sensor MLSL276

[Software de aquisição de dados original do fabricante MLSL276](https://www.wenglor.com/index.php?L=0&id=1148&tx_wsshoploader_pi1[url]=catalog/productDetail.jsf;jsessionid::2fRby-Npc6Ls2NlJg7M-umOeFBqNYwEu6CoA_SAPOSQonLixnzH0JzeTUavcw46w;saplb_*::(J2EE2811920)2811950;;wec-appid::Shop_1000_EXT_EN;;itemKey::MLSL276;;wec-locale::en_US;;ifr::y), navegue para a aba de downloads e selecione SDK_Windows_Linux_weCat3D_1.2.0.zip


[Software de aquisição de dados modificado para o projeto VIMTILO](https://www.dropbox.com/s/vs0yw5p3ww2lnqf/weCat3D_SDK_Windows_QT_C%2B%2B_V_2_1_3.rar?dl=0)

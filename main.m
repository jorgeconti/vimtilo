%% inicializacao
% Na inicialização é possível alterar algumas variáveis que devem ser
% ajustadas para cada ponte.
clc; clearvars; close all;

%% salvar arquivos
salvarPNG = false;
salvarFIG = false;
salvarEPS = false;
salvarpontoM = false;

%% plot de graficos
plotsDesnecessarios = false;
plotreduzido = false;

%% nome do arquivo
nome='ponte_inda_posicao_normal_1mm_laser_alto.mat';
destinoImagens="D:\data\imagens\";

%% Carregar dados da ponte
load(nome)

%% especificações do novo dormente
dormenteComprimento=2700;
dormenteLargura=200;
dormenteDistanciadeCentro=380;
distanciaLongarinas=143;

dormenteMesmaPosicao=1;

%% Configuração do encoder
pulsos_encoder=600;
diametro_eixo_encoder = 12.9398; %12.73;%*0.98379;%mm
dir_encoder=-1;

%% configurações de ajuste

minArea=10; %% area minima do blob para ser considerado nos calculos

threshold=36;%distancia minima de cada centroide para longarinas diagonais dos outros centroides

thresholdColuna=5;%distancia minima em  pixeis  para considerar outros valores como coluna

thresholddiag=thresholdColuna;% remove a colunas verticais deixando apenas as diagonais

pecentual_Largura_Contraventamentos = 2;
%% ajuste de pre filtagem da nuvem de pontos

xoff=1510;
margem_Superior_Longarina_Z=1390;
margem_Inferior_Longarina_Z=1355;
margem_Superior_Via_Z=1200;
margem_Inferior_Via_Z=924;
margem_Superior_Contraventamento_Z=1650;
margem_Inferior_Contraventamento_Z=1450;
margem_Inferior_Z_trilho=900;
margem_Superior_Y=15000;
margem_Inferior_Y=13500;
margem_Superior_X=200;
margem_Inferior_X=-200;
margem_Superior_I=500;
margem_Inferior_I=100;
margem_Inferior_I_trilho=80;


%% variaveis de processo

global maiorY menorY novoMinY novoMaxY maiorX menorX novoMinX novoMaxX x y
x=1;
y=2;
%% espelhamento axial
EX=-1;
EY=-1;
EZ=1;
EI=-1;

%% resolução da nuvem de pontos
reducao=10;
resolucao

%% Gerar mapa de caracteristicas da nuvem de pontos
Y=aux1E*(pi*200)/20000;

X11=aux1X(:,:)+xoff;
Y11=Y(:,:)+10;
Z11=aux1Z(:,:)-7;
I11=aux1I(:,:);

Z11=Z11.*(X11~=xoff);
Y11=Y11.*(X11~=xoff);
I11=I11.*(X11~=xoff);
X11=X11.*(X11~=xoff);

Zsc=Z11<margem_Superior_Contraventamento_Z;
Zic=Z11>margem_Inferior_Contraventamento_Z;
Zsl=Z11<margem_Superior_Longarina_Z;
Zil=Z11>margem_Inferior_Longarina_Z;
Zsv=Z11<margem_Superior_Via_Z;
Ziv=Z11>margem_Inferior_Via_Z;
Zt=Z11>margem_Inferior_Z_trilho;
Ys=Y11<margem_Superior_Y;
Yi=Y11>margem_Inferior_Y;
Xs=X11<margem_Superior_X;
Xi=X11>margem_Inferior_X;
Is=I11<margem_Superior_I;
Ii=I11>margem_Inferior_I;
It=I11>margem_Inferior_I_trilho;

Ymap=Yi&Ys;
Xmap=Xi&Xs;
Zmapc=Zic&Zsc;
Zmapl=Zil&Zsl;
Zmapv=Ziv&Zsv;
ZmapTrilho=Zt&Zsc;
Imap=Ii&Is;
ImapTrilho=It&Is;

dotcloudContraventamento1=Zmapc;
dotcloudLongarina1       =Zmapl;
dotcloudVia1             =Zmapv;
dotcloudcutContraventamento1=Xmap&Zmapc&Imap;
dotcloudcutLongarina1       =Xmap&Zmapl&Imap;
dotcloudcutVia1             =Xmap&Zmapv&Imap;
dotcloudcutTrilho1=Xmap&ZmapTrilho&ImapTrilho;

I22=aux2I(:,:);
Z22=aux2Z(:,:);
X22=aux2X(:,:);
Y22=Y(:,:);
Z22=Z22.*(X22~=0);
Y22=Y22.*(X22~=0);
I22=I22.*(X22~=0);
X22=X22.*(X22~=0);
Zsc=Z22<margem_Superior_Contraventamento_Z;
Zic=Z22>margem_Inferior_Contraventamento_Z;
Zsl=Z22<margem_Superior_Longarina_Z;
Zil=Z22>margem_Inferior_Longarina_Z;
Zsv=Z22<margem_Superior_Via_Z;
Ziv=Z22>margem_Inferior_Via_Z;
Zt=Z22>margem_Inferior_Z_trilho;
Ys=Y22<margem_Superior_Y;
Yi=Y22>margem_Inferior_Y;
Xs=X22<margem_Superior_X;
Xi=X22>margem_Inferior_X;
Is=I22<margem_Superior_I;
Ii=I22>margem_Inferior_I;
It=I22>margem_Inferior_I_trilho;
Ymap=Yi&Ys;
Xmap=Xi&Xs;
Zmapc=Zic&Zsc;
Zmapl=Zil&Zsl;
Zmapv=Ziv&Zsv;
ZmapTrilho=Zt&Zsc;
Imap=Ii&Is;
ImapTrilho=It&Is;

dotcloudContraventamento2=Zmapc;
dotcloudLongarina2       =Zmapl;
dotcloudVia2             =Zmapv;
dotcloudcutContraventamento2=Xmap&Zmapc&Imap;
dotcloudcutLongarina2       =Xmap&Zmapl&Imap;
dotcloudcutVia2             =Xmap&Zmapv&Imap;
dotcloudcutTrilho2=Xmap&ZmapTrilho&ImapTrilho;

clearvars aux1E aux1I aux1T aux1X aux1Z aux2E aux2I aux2T aux2X aux2Z;
clearvars Ymap Yi Ys Xmap Xi Xs Zmapc Zic Zsc Zmapl Zil Zsl Zmapv Ziv Zsv;
clearvars ZmapTrilho Zt Zsc Imap Ii Is ImapTrilho It Is resolusao

%%%%%%%%%%%%%%%%%%%%%%%% vias
vias = [dotcloudVia1 dotcloudVia2];
if(plotsDesnecessarios)
    figure;
    imshow(vias);
end
if(salvarpontoM)
    save(destinoImagens+'vias.m','vias');
end
clearvars dotcloudVia1 dotcloudVia2


%%%%%%%%%%%%%%%%%%%%%%%% longarinas
longarinas = [dotcloudLongarina1 dotcloudLongarina2];
if(plotsDesnecessarios)
    imshow(longarinas);
end
if(salvarpontoM)
    save(destinoImagens+'longarinas.m','longarinas');
end
clearvars dotcloudLongarina1 dotcloudLongarina2


%%%%%%%%%%%%%%%%%%%%%%%% contraventamentos
contraventamentos = [dotcloudContraventamento1 dotcloudContraventamento2];
if(plotsDesnecessarios)
    figure;
    imshow(contraventamentos);
end
if(salvarpontoM)
    save(destinoImagens+'contraventamentos.m','contraventamentos');
end
clearvars dotcloudContraventamento1 dotcloudContraventamento2

XX=[X11 X22];
YY=[Y11 Y22];
ZZ=[Z11 Z22];
II=[I11 I22];

clearvars X11 Y11 Z11 I11 X22 Y22 Z22 I22

dotcloudcut = longarinas | contraventamentos;

maiorX=size(dotcloudcut,2);
menorX=1;
maiorY=size(dotcloudcut,1);
menorY=1;
novoMinY=min(min(Y));
novoMaxY=max(max(Y));

X_dot = XX(dotcloudcut);
Y_dot = YY(dotcloudcut);
Z_dot = ZZ(dotcloudcut);
I_dot = II(dotcloudcut);
X_dot = X_dot(:);
Y_dot = Y_dot(:);
Z_dot = Z_dot(:);
I_dot = I_dot(:);

X_dot_trilho = XX(vias);
Y_dot_trilho = YY(vias);
Z_dot_trilho = ZZ(vias);
I_dot_trilho = II(vias);
X_dot_trilho = X_dot_trilho(:);
Y_dot_trilho = Y_dot_trilho(:);
Z_dot_trilho = Z_dot_trilho(:);
I_dot_trilho = I_dot_trilho(:);

clear Y X Z I;
%% redimencionarImagem2NuvemLimites
% responsável por definir os limites do valor X para determinada altura Z na
% qual todos os cálculos foram realizados, sendo possível assim retornar
% toda a informação obtida na imagem binaria para a nuvem de pontos, e
% assim representando os dados obtidos em valores métricos reais

mlslRangeZ=[300 1500];
mlslWidthX=[250 1350];
mlslAngle=((mlslWidthX(2)/2-mlslWidthX(1)/2)/(mlslRangeZ(2)-mlslRangeZ(1)));
WidthoX=((mlslAngle)*(1680-mlslRangeZ(1)))+((mlslWidthX(1)/2));
WidthoX=WidthoX*2;

offsetReal=-760;
novoMaxX=WidthoX+offsetReal;
novoMinX=-WidthoX+offsetReal;
%%

%% plotarScatterInicial
% realiza o plot tridimensional inicial da nuvem de pontos com
% predefinições implementadas em (inicializacao)

if(plotsDesnecessarios)
    figure('Name','Nuvem de pontos','NumberTitle','off');
    scatter3(EX*X_dot_trilho,EZ*Z_dot_trilho,EY*Y_dot_trilho,1,EI*I_dot_trilho)
    title({'Nuvem de pontos formada com os dados brutos do sensor MLSL276'})
    
    zlabel('Y  [mm]');
    xlabel('X  [mm]');
    ylabel('Z  [mm]');
end

%%% salvar figura
nomeFigura="Nuvem de pontos";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end

%% plotarScatterInicialpos longarinas
% realiza o plot tridimensional da longarina

if(plotsDesnecessarios)
    figure('Name','Nuvem de pontos','NumberTitle','off');
    scatter3(EX*XX(longarinas),EY*YY(longarinas),-EZ*ZZ(longarinas),1,EZ*ZZ(longarinas))
    title({'Parte estrutural da longarina pós filtragem'})
    zlabel('Z  [mm]');
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
end
%%% salvar figura
nomeFigura="Nuvem de pontos longarina filtrada";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end

%% plotarScatterInicialpos filtragem
% realiza o plot tridimensional inicial da nuvem de pontos com
% predefinições implementadas em (inicializacao)

if(plotsDesnecessarios)
    figure('Name','Nuvem de pontos','NumberTitle','off');
    scatter3(EX*X_dot,EZ*Z_dot,EY*Y_dot,1,EI*I_dot)
    c = colorbar;
    c.Label.String = 'Intensidade luminosa';
    title('normal')
    xlabel('X')
    ylabel('Z')
    zlabel('Y')
    title({'Parte estrutural da ponte pós filtragem'})
    
    zlabel('Y  [mm]');
    xlabel('X  [mm]');
    ylabel('Z  [mm]');
end
%%% salvar figura
nomeFigura="Nuvem de pontos filtrada";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end
%% inicializarHblod
% Trabalhando com os blobs gerados

hblob = vision.BlobAnalysis;
hblob.MinimumBlobArea = minArea;

hblob.AreaOutputPort = true;
hblob.CentroidOutputPort = true;
hblob.BoundingBoxOutputPort = true;
hblob.MajorAxisLengthOutputPort = true;
hblob.MinorAxisLengthOutputPort = true;
hblob.OrientationOutputPort = true;
hblob.EccentricityOutputPort = true;
hblob.EquivalentDiameterSquaredOutputPort = true;
hblob.ExtentOutputPort = true;
hblob.PerimeterOutputPort = true;
hblob.MaximumCount = 1000;


%% longarinasdeSuporte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, centroid, bbVertival, maiorEixo, ~, ~] = step(hblob, longarinas);

margem_Superior_Y=max(centroid(:,2));
margem_Inferior_Y=min(centroid(:,2));

%%%coluna 1
c1=round(centroid(:,1)); %encontra os valores mais repetidos q sera considerado uma coluna
b=zeros(size(c1));
for i=1:size(c1,1)
    b(i)=sum(c1>(c1(i)-std(c1)/3)&c1<(c1(i)+std(c1)/3));
end

valorColuna1=round(mean(c1(b(1:round(size(b,1)/2))==mode(b(1:round(size(b,1)/2))))));
valorColuna2=round(mean(c1(logical([zeros(round(size(b,1)/2)-1,1); b(round(size(b,1)/2):end)==mode(b(round(size(b,1)/2):end))]))));

if(valorColuna1<valorColuna2) %determina o vetor logico sempre mantendo a coluna 1 como a longarina da esquerda
    vecColuna1=centroid<(valorColuna1+thresholdColuna)&centroid>(valorColuna1-thresholdColuna);
    vecColuna2=centroid<(valorColuna2+thresholdColuna)&centroid>(valorColuna2-thresholdColuna);
else
    aux=valorColuna1;
    valorColuna1=valorColuna2;
    valorColuna2=aux;
    vecColuna1=centroid<(valorColuna1+thresholdColuna)&centroid>(valorColuna1-thresholdColuna);
    vecColuna2=centroid<(valorColuna2+thresholdColuna)&centroid>(valorColuna2-thresholdColuna);
end

vecColuna1(:,2)=0;
vecColuna2(:,2)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EixovcE=double(maiorEixo).*vecColuna1(:,1);
mediaEixovcE=mean(EixovcE(EixovcE~=0));
rangeEixoE=max(EixovcE)-min(EixovcE(EixovcE~=0));
limiteSuperiorEixoE=rangeEixoE/4;
vecColunaEixoE=(EixovcE-mediaEixovcE)<limiteSuperiorEixoE&(EixovcE-mediaEixovcE)>-mediaEixovcE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EixovcD=double(maiorEixo).*vecColuna2(:,1);
mediaEixovcD=mean(EixovcD(EixovcD~=0));
rangeEixoD=max(EixovcD)-min(EixovcD(EixovcD~=0));
limiteSuperiorEixoD=rangeEixoD/4;
vecColunaEixoD=(EixovcD-mediaEixovcD)<limiteSuperiorEixoD&(EixovcD-mediaEixovcD)>-mediaEixovcD;

xc = centroid(:, 1).*(vecColunaEixoE|vecColunaEixoD);
yc = centroid(:, 2).*(vecColunaEixoE|vecColunaEixoD);


xc1=xc(xc~=0&vecColunaEixoE);
yc1=yc(yc~=0&vecColunaEixoE);
coluna1 = fit( yc1,xc1 ,  'poly1' );% reta gerada para longarina esquerda

xc2=xc(xc~=0&vecColunaEixoD);
yc2=yc(yc~=0&vecColunaEixoD);
coluna2 = fit( yc2 ,xc2,  'poly1' );% reta encontrada para longarina da direita

centroidSortY=sortrows(centroid(vecColuna1(:,1)|vecColuna2(:,1),:),2);

ylateral=ones(1,size(dotcloudcut, 1));
ylateral=ylateral.*(1:size(dotcloudcut, 1));
xlateralEsq = coluna1.p1*ylateral + coluna1.p2;
xlateralDir = coluna2.p1*ylateral + coluna2.p2;
xlateralEsqInf=xlateralEsq-thresholdColuna;  % aqui é considerado uma margem de largura para a longarina
xlateralEsqSup=xlateralEsq+thresholdColuna;
xlateralDirInf=xlateralDir-thresholdColuna;
xlateralDirSup=xlateralDir+thresholdColuna;
if(plotsDesnecessarios)
    figure('Name','Imagem binaria','NumberTitle','off');
    imshow(dotcloudcut)
end
%%% salvar figura
nomeFigura="Imagem binaria";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end
if(plotsDesnecessarios)
    figure('Name','Imagem binaria com destaque estrutural','NumberTitle','off');
    imshow(dotcloudcut)
    hold on
    % plot(xlateral,ylateral,'r')
    p=plot(xlateralEsqInf,ylateral,'-.b',xlateralEsqSup,ylateral,'-.b',xlateralDirInf,ylateral,'-.b',xlateralDirSup,ylateral,'-.b');
    p(1).LineWidth = 1.3;
    p(2).LineWidth = 1.3;
    p(3).LineWidth = 1.3;
    p(4).LineWidth = 1.3;
    
    p=plot(xc1,yc1,'b*');
    p.LineWidth = 1.3;
    p=plot(xc2,yc2,'b*');
    p.LineWidth = 1.3;
    title('normal')
    title({'Realce dos pontos utilizados para encontrar as Longarinas de sustentação'})
    legend('Longarinas de sustentação')
    hold off
end

%%% salvar figura
nomeFigura="Imagem binaria com destaque estrutural";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end
%% redefinirDotcloudApenasDiagonais
% dotclouddiag=zeros(size(dotcloudcut));%dotclouddiag é gerada para realizar o tratamento das retas diagonais
% que são partes estruturais da ponte que não devem manter contato com os
% dormentes
%
% assim na criação desse mapa de diagonais e removida as longarinas
% encontradas anteriormente
% dotclouddiag(:,(valorColuna1+thresholddiag):(valorColuna2-thresholddiag))= dotcloudcut(:,(valorColuna1+thresholddiag):(valorColuna2-thresholddiag));
% dotclouddiag=dotclouddiag>0.5; % apenas para torna a variável como logical
%
% parte complexa que remove as longarinas laterais da parte central que
% possui as estruturas diagonais

dotclouddiag=zeros(size(dotcloudcut));
dotclouddiag=logical(dotclouddiag);
for i=1:size(dotcloudcut,1)
    for j=min(ceil(xlateralEsqSup)):max(floor(xlateralDirInf))
        dotclouddiag(i,j)=dotcloudcut(i,j)*(xlateralEsqSup(i)<j&xlateralDirInf(i)>j);
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tratarLongarinasDiagonais
% responsável por trabalhar com os centroides das longarinas diagonais e
% trata-los definindo quais são unidos e determinando a melhor reta para cada
% centroide, ao final gerando retas com offset ao redor do eixo definido
% para representas a espessura das longarinas diagonais
% ajustar cada linha diagonal para iniciar e finalizar junto com as


%refazendo o hblod agora para o mapa apenas exclusivamente com diagonais
[area, centroide, bb, maiorEixo, menorEixo, angulo] = step(hblob, contraventamentos);

entreLongarinas = centroide(:,1)>(valorColuna1+thresholddiag/2) & centroide(:,1)<(valorColuna2-thresholddiag/2);
area=area(entreLongarinas);
centroide=centroide(entreLongarinas,:);
bb=bb(entreLongarinas,:);
maiorEixo=maiorEixo(entreLongarinas);
menorEixo=menorEixo(entreLongarinas);
angulo=angulo(entreLongarinas);

limiteMinimodeArea=area>mean(area)*0.05; % remove blobs com area muito pequenas
area=area(limiteMinimodeArea);
centroide=centroide(limiteMinimodeArea,:);
bb=bb(limiteMinimodeArea,:);
maiorEixo=maiorEixo(limiteMinimodeArea);
menorEixo=menorEixo(limiteMinimodeArea);
angulo=angulo(limiteMinimodeArea);


xblob=ones(size(centroide,1),size(dotclouddiag, 2));
xblob=xblob.*(1:size(dotclouddiag, 2));

xc = centroide(:, 1);
yc = centroide(:, 2);
m = -tan(angulo); % Esse sinal de - aparece pois o referencial da imagem é com o Y para baixo

yblob = m.*xblob - m.*xc + yc;

%% determinar matriz retas com a informacao de onde estao localizadas todas

dim=size(xblob,1);
retas=zeros(dim,dim);
for cont = 1:dim
    d=min(sqrt((xblob-xc(cont,1)).^2 + (yblob-yc(cont,1)).^2),[],2);
    retas(:,cont)=d<threshold;
end

%limpa falso eclipse no mesmo eixo
check=retas==retas';
retas=retas.*check;
redundancia = ones(1,size(retas,2)) - sum(tril(retas,-1),1); %elimina valores duplicados de retas
retas=retas(:,redundancia==1);%remove colunas da matriz q sao redundantes para polpar trabalho na hora de plotar
%%

angulo=angulo(redundancia'==1);% remove os angulos q nao seram mais utilizados

nblobs = length(angulo);
xblob = 1:size(dotclouddiag, 2); %Para desenhar a reta proxima à inclinação da longarina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

espessura_da_linha=2.5;
figure('Name','Imagem binaria','NumberTitle','off');

imshow(dotcloudcut);
title({'Imagem binaria com longarinas e dormentação'})
margem=10;
hold on
clear m yblob;
m=zeros(1,nblobs);

condicaoFimReta=xblob>(valorColuna1-thresholdColuna)&xblob<(valorColuna2+thresholdColuna);
maiorYY=0;
menorYY=size(dotcloudcut,1);
contDiagAclive=1;
contDiagDeclive=1;

indexEsq=xlateralEsq>0&xlateralEsq<size(dotcloudcut,2);
indeDir=xlateralDir>0&xlateralDir<size(dotcloudcut,2);

for cont = 1:nblobs
    xc = centroide(retas(:,cont)==1, 1);
    yc = centroide(retas(:,cont)==1, 2);
    
    if(size(xc,1)>1)
        reta1 = fit( xc,yc ,  'poly1' );
        m(cont)=reta1.p1;
    else
        m(cont) = 0;
    end
end
media_angulo_pos=mean(m(m>0));
media_angulo_neg=mean(m(m<0));


for cont = 1:nblobs
    xc = centroide(retas(:,cont)==1, 1);
    yc = centroide(retas(:,cont)==1, 2);
    if(size(xc,1)==1)
        if(-tan(angulo(cont))>0)
            m(cont)=media_angulo_pos;
        else
            m(cont)=media_angulo_neg;
        end
    end
    
    if size(xc,1)>1
        reta2 = fit( xc,yc ,  'poly1' ); % reta gerada para longarina esquerda
        yblob = reta2.p1*xblob + reta2.p2;
    else
        continue %% nenhuma reta com apenas um ponto sera traçada
        yblob =  m(cont)*xblob - m(cont)*xc + yc; %Para desenhar a reta proxima à inclinação da longarina
    end
    
    if(maiorYY<max(max(yblob(:,condicaoFimReta))))
        maiorYY=max(max(yblob(:,condicaoFimReta)));
    end
    if(menorYY>min(min(yblob(:,condicaoFimReta))))
        menorYY=min(min(yblob(:,condicaoFimReta)));
    end
    yblob_menor=mean(yblob(:,:),1)-ceil(mean(menorEixo))*pecentual_Largura_Contraventamentos;
    yblob_maior=mean(yblob(:,:),1)+ceil(mean(menorEixo))*pecentual_Largura_Contraventamentos;
    [~,peme1,peme2,~] =  proximidade(xlateralEsqSup(indexEsq),ylateral(indexEsq), xblob, yblob_menor);
    [~,pdme1,pdme2,~] =  proximidade(xlateralDirInf(indeDir),ylateral(indeDir), xblob, yblob_menor);
    
    [~,pema1,pema2,~] =  proximidade(xlateralEsqSup(indexEsq),ylateral(indexEsq), xblob, yblob_maior);
    [~,pdma1,pdma2,~] =  proximidade(xlateralDirInf(indeDir),ylateral(indeDir), xblob, yblob_maior);
    
    limite_diagonais_menor=xblob>xlateralEsqSup(peme1)&xblob<xlateralDirInf(pdme1);
    limite_diagonais_maior=xblob>xlateralEsqSup(pema1)&xblob<xlateralDirInf(pdma1);
    %     limite_diagonais=xblob>Xeigual&xblob<Xdigual;
    %         plot(xlateralEsqSup(peme1),ylateral(peme1), '*g');
    %         plot(xlateralDirInf(pdme1),ylateral(pdme1), '*g');
    %         plot(xlateralEsqSup(pema1),ylateral(pema1), '*r');
    %         plot(xlateralDirInf(pdma1),ylateral(pdma1), '*r');
    plot(xblob(peme2),yblob_menor(peme2), '*g');
    plot(xblob(pdme2),yblob_menor(pdme2), '*g');
    plot(xblob(pema2),yblob_maior(pema2), '*r');
    plot(xblob(pdma2),yblob_maior(pdma2), '*r');
    if( m(cont)>0)
        p=plot(xc, yc, '*g');
        p.LineWidth = espessura_da_linha;
        p=plot(xblob(limite_diagonais_menor), yblob_menor(limite_diagonais_menor), 'g');
        p.LineWidth = espessura_da_linha;
        p=plot(xblob(limite_diagonais_maior), yblob_maior(limite_diagonais_maior), 'g');
        p.LineWidth = espessura_da_linha;
        Longarina_Diagonal(:,1:sum(limite_diagonais_menor),contDiagDeclive,1,1)=[ xblob(limite_diagonais_menor); yblob_menor(limite_diagonais_menor)];
        Longarina_Diagonal(:,1:sum(limite_diagonais_maior),contDiagDeclive,2,1)=[ xblob(limite_diagonais_maior); yblob_maior(limite_diagonais_maior)];
        contDiagDeclive=contDiagDeclive+1;
    end
    if(m(cont)<0)
        
        p=plot(xc, yc, '*r');
        p.LineWidth = espessura_da_linha;
        p=plot(xblob(limite_diagonais_menor), yblob_menor(limite_diagonais_menor), '--r');
        p.LineWidth = espessura_da_linha;
        p=plot(xblob(limite_diagonais_maior), yblob_maior(limite_diagonais_maior), '--r');
        p.LineWidth = espessura_da_linha;
        Longarina_Diagonal(:,1:sum(limite_diagonais_menor),contDiagAclive,1,2)= [  xblob(limite_diagonais_menor); yblob_menor(limite_diagonais_menor)];
        Longarina_Diagonal(:,1:sum(limite_diagonais_maior),contDiagAclive,2,2)= [  xblob(limite_diagonais_maior); yblob_maior(limite_diagonais_maior)];
        contDiagAclive=contDiagAclive+1;
    end
    
    % plot(xblob, yblob, 'w');
end

%% plota longarinas verticais
% condicaoFimLongarinaEsq=xlateralEsq>menorYY*0.1&xlateralEsq<maiorYY*1.1;
% condicaoFimLongarinaDir=xlateralDir>menorYY*0.1&xlateralDir<maiorYY*1.1;
condicaoFimLongarinaEsq=xlateralEsq~=0;
condicaoFimLongarinaDir=xlateralDir~=0;

xcs = centroid(vecColuna1|vecColuna2, 1);
ycs = centroid(vecColuna1|vecColuna2, 2);
p=plot( xcs,ycs,'*b');
p.LineWidth = espessura_da_linha;

p=plot( xlateralEsqInf(condicaoFimLongarinaEsq),ylateral(condicaoFimLongarinaEsq),'-.b');
p.LineWidth = espessura_da_linha;
p=plot( xlateralEsqSup(condicaoFimLongarinaEsq),ylateral(condicaoFimLongarinaEsq),'-.b');
p.LineWidth = espessura_da_linha;
p=plot( xlateralDirInf(condicaoFimLongarinaDir),ylateral(condicaoFimLongarinaDir),'-.b');
p.LineWidth = espessura_da_linha;
p=plot( xlateralDirSup(condicaoFimLongarinaDir),ylateral(condicaoFimLongarinaDir),'-.b');
p.LineWidth = espessura_da_linha;

Longarina_Vertical(:,1:size(xlateralEsqInf(condicaoFimLongarinaEsq),2),1)=[ xlateralEsqInf(condicaoFimLongarinaEsq);ylateral(condicaoFimLongarinaEsq)];
Longarina_Vertical(:,1:size(xlateralEsqSup(condicaoFimLongarinaEsq),2),2)=[ xlateralEsqSup(condicaoFimLongarinaEsq);ylateral(condicaoFimLongarinaEsq)];
Longarina_Vertical(:,1:size(xlateralDirInf(condicaoFimLongarinaDir),2),3)=[ xlateralDirInf(condicaoFimLongarinaDir);ylateral(condicaoFimLongarinaDir)];
Longarina_Vertical(:,1:size(xlateralDirSup(condicaoFimLongarinaDir),2),4)=[ xlateralDirSup(condicaoFimLongarinaDir);ylateral(condicaoFimLongarinaDir)];
%%

dormenteComprimentoPixel=menorX+((dormenteComprimento- menorX)/( novoMaxX-novoMinX))*(maiorX- menorX);
dormenteLarguraPixel=round(menorY+((dormenteLargura- menorY)/( novoMaxY-novoMinY))*(maiorY-menorY));
dormenteDistanciadeCentroPixel=round(menorY+((dormenteDistanciadeCentro- menorY)/( novoMaxY-novoMinY))*(maiorY- menorY));
c2=round(centroidSortY(:,2)); %encontra os valores mais repetidos q sera considerado uma coluna
b2=[zeros(size(c2)) zeros(size(c2))];
for i=1:size(c2,1)
    b2(i,1)=sum(c2>(c2(i)-std(c2)/(4*reducao))&c2<(c2(i)+std(c2)/(4*reducao)));
    b2(i,2)=i;
end


[Mei,inicio,~,~] =  proximidade(ylateral,0,min(centroidSortY(:,2)),0);
[Mef,fim,~,~] =  proximidade(ylateral,0,max(centroidSortY(:,2)),0);

dormenteCount=1;
dormente=zeros(sum(b2(:,1)==2)/2+sum(b2(:,1)==1),4);
pula=0;

comprimento_dormentacao=abs(margem_Superior_Y)-abs(margem_Inferior_Y);
quantidade_dormentes=round(comprimento_dormentacao/dormenteDistanciadeCentro)+1;

contDormente=1;
for i=1:size(c2,1)
    
    if(pula==1)
        pula=0;
    else
        if(b2(i,1)==2)
            
            %
            if(dormenteMesmaPosicao~=0)
                dormenteAngulo=(centroidSortY(i,y)-centroidSortY(i+1,y))/(centroidSortY(i,x)-centroidSortY(i+1,x))  ;
                passo=round((centroidSortY(i,y)+centroidSortY(i+1,y))/2)+dormenteDistanciadeCentroPixel-round(dormenteLarguraPixel);
                dormenteCentroid=[(xlateralEsq(passo)+(xlateralDir(passo)-xlateralEsq(inicio))/2) ylateral(passo) ];
                dormente(dormenteCount,:)=[dormenteCentroid(1)-dormenteComprimentoPixel/2 dormenteCentroid(1)+dormenteComprimentoPixel/2 dormenteCentroid(2)-dormenteLarguraPixel/2 dormenteCentroid(2)+dormenteLarguraPixel/2];
                limitesX=xblob>dormente(dormenteCount,1)&xblob<dormente(dormenteCount,2);
                yblobDS =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,3));
                yblobDI =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,4));
                %             plot(xblob(limitesX),yblobDS(limitesX),'y')
                %             plot(xblob(limitesX),yblobDI(limitesX),'y')
                dormente_salvo(:,1:sum(limitesX),contDormente,1)=[xblob(limitesX);yblobDS(limitesX)];
                dormente_salvo(:,1:sum(limitesX),contDormente,2)=[xblob(limitesX);yblobDI(limitesX)];
                if(dormenteAngulo<0)
                    dormente_pontos(:,:,contDormente,3)=[min(xblob(limitesX)) max(xblob(limitesX)) min(xblob(limitesX)) max(xblob(limitesX)) ; max(yblobDS(limitesX)) min(yblobDS(limitesX))  max(yblobDI(limitesX)) min(yblobDI(limitesX))];
                else
                    dormente_pontos(:,:,contDormente,3)=[min(xblob(limitesX)) max(xblob(limitesX)) min(xblob(limitesX)) max(xblob(limitesX)) ; min(yblobDS(limitesX)) max(yblobDS(limitesX))  min(yblobDI(limitesX)) max(yblobDI(limitesX))];
                end
                contDormente=contDormente+1;
                pula=1;
            end
        end
        if(b2(i,1)==1)
            if(dormenteCount>1)
                if(dormenteMesmaPosicao~=0)
                    dormenteCentroid=[dormenteCentroid(1) (centroidSortY(i,2)+dormenteDistanciadeCentroPixel-dormenteLarguraPixel)];
                    dormente(dormenteCount,:)=[dormenteCentroid(1)-dormenteComprimentoPixel/2 dormenteCentroid(1)+dormenteComprimentoPixel/2 dormenteCentroid(2)-dormenteLarguraPixel/2 dormenteCentroid(2)+dormenteLarguraPixel/2];
                    limitesX=xblob>dormente(dormenteCount,1)&xblob<dormente(dormenteCount,2);
                    yblobDS =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,3));
                    yblobDI =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,4));
                    %                 plot(xblob(limitesX),yblobDS(limitesX),'y')
                    %                 plot(xblob(limitesX),yblobDI(limitesX),'y')
                end
            end
        end
        dormenteCount = dormenteCount + 1;
    end
end
% end

if(dormenteMesmaPosicao==0)
    for i=1:quantidade_dormentes
        dormenteAngulo=EY*mean([coluna1.p1 coluna2.p1]) ;
        %     passo=inicio+dormenteDistanciadeCentroPixel*(i-1);%i-dormenteLarguraPixel;
        passo=inicio+dormenteDistanciadeCentroPixel*i-round(dormenteLarguraPixel*1.4);
        dormenteCentroid=[(xlateralEsq(passo)+(xlateralDir(passo)-xlateralEsq(inicio))/2) ylateral(passo) ];
        dormente(dormenteCount,:)=[dormenteCentroid(1)-dormenteComprimentoPixel/2 dormenteCentroid(1)+dormenteComprimentoPixel/2 dormenteCentroid(2)-dormenteLarguraPixel/2 dormenteCentroid(2)+dormenteLarguraPixel/2];
        limitesX=xblob>dormente(dormenteCount,1)&xblob<dormente(dormenteCount,2);
        yblobDS =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,3));
        yblobDI =  (dormenteAngulo*xblob - dormenteAngulo*dormenteCentroid(1) + dormente(dormenteCount,4));
        %             plot(xblob(limitesX),yblobDS(limitesX),'y')
        %             plot(xblob(limitesX),yblobDI(limitesX),'y')
        dormente_salvo(:,1:sum(limitesX),i,1)=[xblob(limitesX);yblobDS(limitesX)];
        dormente_salvo(:,1:sum(limitesX),i,2)=[xblob(limitesX);yblobDI(limitesX)];
        if(dormenteAngulo<0)
            dormente_pontos(:,:,i,3)=[min(xblob(limitesX)) max(xblob(limitesX)) min(xblob(limitesX)) max(xblob(limitesX)) ; max(yblobDS(limitesX)) min(yblobDS(limitesX))  max(yblobDI(limitesX)) min(yblobDI(limitesX))];
        else
            dormente_pontos(:,:,i,3)=[min(xblob(limitesX)) max(xblob(limitesX)) min(xblob(limitesX)) max(xblob(limitesX)) ; min(yblobDS(limitesX)) max(yblobDS(limitesX))  min(yblobDI(limitesX)) max(yblobDI(limitesX))];
        end
    end
end
%% contorno dos dormentes
for di=1:size(dormente_pontos,3)
    han=[[dormente_pontos(x,1,di,3);dormente_pontos(x,2,di,3);dormente_pontos(x,4,di,3);dormente_pontos(x,3,di,3);dormente_pontos(x,1,di,3)] [dormente_pontos(y,1,di,3);dormente_pontos(y,2,di,3);dormente_pontos(y,4,di,3);dormente_pontos(y,3,di,3);dormente_pontos(y,1,di,3)]];
    p=plot(han(:,1),han(:,2),':y');
    p.LineWidth = espessura_da_linha;
    clear han
end

%%% salvar figura
nomeFigura="Imagem binaria com longarinas e dormentação";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end


figure('Name','Dormentação','NumberTitle','off');
title({'Dormetação com entalhes '})
hold on

CDP=ones(size(dormente_salvo,3),1);
%% identifica e armazena interseções entre dormentes e longarinas diagonais
for dq=1:2
    for di=1:size(dormente_salvo,3)
        %plot(dormente_salvo(x,:,di,dq), dormente_salvo(y,:,di,dq),'k')
        
        for ldq=1:2
            for ldl=1:2
                for ld=1:size(Longarina_Diagonal,3)
                    pld(:,ldl,ldq)=Longarina_Diagonal(y,:,ld,ldl,ldq)~=0;
                    [Me,pos1,pos2,Xigual] =  proximidade(Longarina_Diagonal(x,pld(:,ldl,ldq),ld,ldl,ldq),Longarina_Diagonal(y,pld(:,ldl,ldq),ld,ldl,ldq), dormente_salvo(x,:,di,dq), dormente_salvo(y,:,di,dq));
                    if(Me<1)
                        %       plot(Longarina_Diagonal(x,pld(:,ldl,ldq),ld,ldl,ldq),Longarina_Diagonal(y,pld(:,ldl,ldq),ld,ldl,ldq),'b')
                        dormente_pontos(:,CDP(di),di,1)=[dormente_salvo(x,pos2,di,dq);dormente_salvo(y,pos2,di,dq)];
                        %        plot(dormente_pontos(x,CDP(di),di,1),dormente_pontos(y,CDP(di),di,1),'r*')
                        Qual_Longarina(:,di,dq,ld, ldl, ldq)=[dormente_salvo(x,pos2,di,dq);dormente_salvo(y,pos2,di,dq)];
                        CDP(di)=CDP(di)+1;
                    end
                end
            end
        end
    end
end
%% identifica e armazena interseções entre dormentes e longarinas verticais
CDP=ones(size(dormente_salvo,3),1);
for dq=1:2
    for di=1:size(dormente_salvo,3)
        %plot(dormente_salvo(x,:,di,dq), dormente_salvo(y,:,di,dq),'k')
        for lv=1:size(Longarina_Vertical,3)
            plv(:,lv)=Longarina_Vertical(y,:,lv)~=0;
            [Me,pos1,pos2,Xigual] =  proximidade(Longarina_Vertical(x,plv(:,lv),lv),Longarina_Vertical(y,plv(:,lv),lv), dormente_salvo(x,:,di,dq), dormente_salvo(y,:,di,dq));
            if(Me<50)
                %                 plot(Longarina_Vertical(x,plv(:,lv),lv),Longarina_Vertical(y,plv(:,lv),lv),'k')
                dormente_pontos(:,CDP(di),di,2)=[dormente_salvo(x,pos2,di,dq);dormente_salvo(y,pos2,di,dq)];
                CDP(di)=CDP(di)+1;
            end
        end
    end
end

clear pld

%% plota os pedaços de longarinas diagonais que passam por dentro do dormente
espessura_da_linha=2.5;
pontoextra=1;
contaux1=0;
contaux2=0;
for dq=1:2
    for di=1:size(Qual_Longarina,2)
        for ldq=1:2
            for ldl=1:2
                for ld=1:size(Longarina_Diagonal,3)
                    pld(:,ldl,ldq)=Longarina_Diagonal(y,:,ld,ldl,ldq)~=0;
                    %% longarina diagonal que atraveça o dormente
                    if(Qual_Longarina(y,di,1,ld, ldl, ldq)~=0 && Qual_Longarina(y,di,2,ld, ldl, ldq)~=0)
                        entre_dormente(:,ld,ldl,ldq)=Longarina_Diagonal(y,:,ld,ldl,ldq)<Qual_Longarina(y,di,2,ld, ldl, ldq)&Longarina_Diagonal(y,:,ld,ldl,ldq)>Qual_Longarina(y,di,1,ld, ldl, ldq);
                        p=plot(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),'--r');
                        big=size(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),2);
                        contaux1=contaux1+1;
                        dormente_entalhe_diagonal(:,1:big,contaux1,1)=[Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq);Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq)];
                        p.LineWidth = espessura_da_linha;
                        
                    end
                    %% longarina diagonal que passa apenas na parte superior do dormente
                    if(Qual_Longarina(y,di,1,ld, ldl, ldq)~=0 && Qual_Longarina(y,di,2,ld, ldl, ldq)==0)
                        entre_dormente(:,ld,ldl,ldq)=Longarina_Diagonal(y,:,ld,ldl,ldq)>Qual_Longarina(y,di,1,ld, ldl, ldq);
                        p=plot(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),'g');
                        big=size(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),2);
                        contaux2=contaux2+1;
                        dormente_entalhe_diagonal(:,1:big,contaux2,2)=[Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq);Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq)];
                        
                        if(sum(entre_dormente(:,ld,ldl,ldq))>0)
                            p.LineWidth = espessura_da_linha;
                            if(ldq==1)
                                dormente_pontos(:,pontoextra,di,4)=[Longarina_Diagonal(x,sum(pld(:,ldl,ldq)),ld,ldl,ldq);Longarina_Diagonal(y,sum(pld(:,ldl,ldq)),ld,ldl,ldq)];
                                pontoextra=pontoextra+1;
                            else
                                dormente_pontos(:,pontoextra,di,4)=[Longarina_Diagonal(x,pld(1,ldl,ldq),ld,ldl,ldq);Longarina_Diagonal(y,pld(1,ldl,ldq),ld,ldl,ldq)];
                                pontoextra=pontoextra+1;
                            end
                        end
                    end
                    %% longarina diagonal que passa apenas na parte superior inferior
                    if(Qual_Longarina(y,di,1,ld, ldl, ldq)==0 && Qual_Longarina(y,di,2,ld, ldl, ldq)~=0)
                        entre_dormente(:,ld,ldl,ldq)=Longarina_Diagonal(y,:,ld,ldl,ldq)<Qual_Longarina(y,di,2,ld, ldl, ldq)&Longarina_Diagonal(y,:,ld,ldl,ldq)~=0;
                        p=plot(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),'g');
                        big=size(Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq),2);
                        contaux2=contaux2+1;
                        dormente_entalhe_diagonal(:,1:big,contaux2,2)=[Longarina_Diagonal(x,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq);Longarina_Diagonal(y,entre_dormente(:,ld,ldl,ldq),ld,ldl,ldq)];
                        
                        if(sum(entre_dormente(:,ld,ldl,ldq))>0)
                            p.LineWidth = espessura_da_linha;
                            if(ldq==1)
                                dormente_pontos(:,pontoextra,di,4)=[Longarina_Diagonal(x,pld(1,ldl,ldq),ld,ldl,ldq);Longarina_Diagonal(y,pld(1,ldl,ldq),ld,ldl,ldq)];
                                pontoextra=pontoextra+1;
                            else
                                dormente_pontos(:,pontoextra,di,4)=[Longarina_Diagonal(x,sum(pld(:,ldl,ldq)),ld,ldl,ldq);Longarina_Diagonal(y,sum(pld(:,ldl,ldq)),ld,ldl,ldq)];
                                pontoextra=pontoextra+1;
                            end
                        end
                    end
                end
            end
        end
    end
end


%% contorno dos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,1,di,3);dormente_pontos(x,2,di,3);dormente_pontos(x,4,di,3);dormente_pontos(x,3,di,3);dormente_pontos(x,1,di,3)] [dormente_pontos(y,1,di,3);dormente_pontos(y,2,di,3);dormente_pontos(y,4,di,3);dormente_pontos(y,3,di,3);dormente_pontos(y,1,di,3)]];
    p=plot(han(:,1),han(:,2),':k');
    p.LineWidth = espessura_da_linha;
    clear han
end
%% realce das longarinas verticais nos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,5,di,2);dormente_pontos(x,1,di,2)] [dormente_pontos(y,5,di,2);dormente_pontos(y,1,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    p.LineWidth = espessura_da_linha;
    clear han
    han=[[dormente_pontos(x,2,di,2);dormente_pontos(x,6,di,2)] [dormente_pontos(y,2,di,2);dormente_pontos(y,6,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    p.LineWidth = espessura_da_linha;
    clear han
end
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,7,di,2);dormente_pontos(x,3,di,2)] [dormente_pontos(y,7,di,2);dormente_pontos(y,3,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    p.LineWidth = espessura_da_linha;
    clear han
    han=[[dormente_pontos(x,4,di,2);dormente_pontos(x,8,di,2)] [dormente_pontos(y,4,di,2);dormente_pontos(y,8,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    p.LineWidth = espessura_da_linha;
    clear han
end
%% plotar todos os pontos importantes

for i=1:size(dormente_pontos,3)
    for j=1:size(dormente_pontos,2)
        if(dormente_pontos(y,j,i,1)~=0)
            plot(dormente_pontos(x,j,i,1),dormente_pontos(y,j,i,1),'r*')
        end
        if(dormente_pontos(y,j,i,2)~=0)
            plot(dormente_pontos(x,j,i,2),dormente_pontos(y,j,i,2),'b*')
        end
        if(dormente_pontos(y,j,i,3)~=0)
            plot(dormente_pontos(x,j,i,3),dormente_pontos(y,j,i,3),'k*')
        end
        if(dormente_pontos(y,j,i,4)~=0)
            plot(dormente_pontos(x,j,i,4),dormente_pontos(y,j,i,4),'g*')
        end
    end
end
xlabel('Pixel');
ylabel('Pixel');
if(EY<0)
    axis ij
else
    axis xy
end

%%% salvar figura
nomeFigura="Dormetação com entalhes ";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end

%% Converter as medidas dos pontos de interseção nos dormentes para tamanho real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%XXXXXXXXXXXx
dormente(:,1:2)=convert_to_real(dormente(:,1:2),x);
dormente_pontos(x,:,:,:)=convert_to_real(dormente_pontos(x,:,:,:),x);
Longarina_Vertical(x,:,:)=convert_to_real(Longarina_Vertical(x,:,:),x);
limitInfXplot=convert_to_real(min(centroid(:,1)),x)-600;
limitSupXplot=convert_to_real(max(centroid(:,1)),x)+400;
dormente_entalhe_diagonal(x,:,:,:)=convert_to_real(dormente_entalhe_diagonal(x,:,:,:),x);
Longarina_Diagonal(x,:,:,:,:)=convert_to_real(Longarina_Diagonal(x,:,:,:,:),x);
centroide(:,x)=convert_to_real(centroide(:,x),x);
%%%%%%%%%%%%YYYYYYYY
valor_condicional=convert_to_real(0,y);
dormente(:,3:4)=convert_to_real(dormente(:,3:4),y)*EY;
dormente_entalhe_diagonal(y,:,:,:)=convert_to_real(dormente_entalhe_diagonal(y,:,:,:),y)*EY;
dormente_entalhe_diagonal(dormente_entalhe_diagonal==valor_condicional*EY)=0;
dormente_pontos(y,:,:,:)=convert_to_real(dormente_pontos(y,:,:,:),y)*EY;
dormente_pontos(dormente_pontos==valor_condicional*EY)=0;
Longarina_Vertical(y,:,:)=convert_to_real(Longarina_Vertical(y,:,:),y)*EY;
Longarina_Vertical(Longarina_Vertical==valor_condicional*EY)=0;
Qual_Longarina(y,:,:,:,:,:)=convert_to_real(Qual_Longarina(y,:,:,:,:,:),y);
Qual_Longarina(Qual_Longarina==valor_condicional)=0;
Longarina_Diagonal(y,:,:,:,:,:)=convert_to_real(Longarina_Diagonal(y,:,:,:,:,:),y);
Longarina_Diagonal(Longarina_Diagonal==valor_condicional)=0;
centroide(:,y)=convert_to_real(centroide(:,y),y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotNuvemdePontosComLongarinasTracadas
% responsável por unir os dados de entrada da nuvem com os calculados em
% imagens binarias e assim assimilar o resultado final realizando a
% comparação em três plots;
% 1- nuvem de pontos com recortes, apresentando a ponte e a linha férrea
% 2- nuvem de pontos similar a 1 sem a linha férrea
% 3- nuvem de pontos similar a 2 com a inclusão de partes estruturais da
% ponte realçadas por linhas definidas em (longarinasdeSuporte) e
% (tratarLongarinasDiagonais)
%%
%%% salvar figura

espessura_da_linha=1.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plotreduzido)
    margem_Inferior_Y = -2.2820e+04
    margem_Superior_Y = -2.6700e+04
else
    margem_Inferior_Y=max(max(EY*Y_dot_trilho));
    margem_Superior_Y=min(min(EY*Y_dot_trilho));
end
figure('Name','Visualização das etapas de filtragem, criação das partes estruturais da ponte e entalhes da dos dormentes','NumberTitle','off');
ax1=subplot(1,3,1);
scatter(EX*X_dot_trilho,EY*Y_dot_trilho,1,EI*Z_dot_trilho)
title({'Nuvem de pontos da ponte';'com os dormentes'})
xlabel('X')
ylabel('Y')
zlabel('Z')
if(dir_encoder>0)
    axis([limitInfXplot limitSupXplot margem_Inferior_Y margem_Superior_Y ])
else
    axis([limitInfXplot limitSupXplot  margem_Superior_Y margem_Inferior_Y])
end
if(EY>0)
    axis ij
else
    axis xy
end
c = colorbar;
c.Label.String = 'Z  [mm]';
xlabel({'X  [mm]';'A'});
ylabel('Y  [mm]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2);
scatter(EX*X_dot,EY*Y_dot ,1,EI*Z_dot)
title({'Nuvem de pontos da ponte';'com estrutura reconstruída'})
xlabel('X')
ylabel('Y')
zlabel('Z')
if(dir_encoder>0)
    axis([limitInfXplot limitSupXplot margem_Inferior_Y margem_Superior_Y ])
else
    axis([limitInfXplot limitSupXplot  margem_Superior_Y margem_Inferior_Y])
end
if(EY>0)
    axis ij
else
    axis xy
end
c = colorbar;
c.Label.String = 'Z  [mm]';
xlabel({'X  [mm]';'B'});
ylabel('Y  [mm]');
hold on
for i=1:size(Longarina_Vertical,3)
    p0=Longarina_Vertical(y,:,i)~=0;
    p=plot( Longarina_Vertical(x,p0,i),Longarina_Vertical(y,p0,i),'-.b');
    p.LineWidth = espessura_da_linha;
end

for i=1:size(Longarina_Diagonal,3)
    p1=Longarina_Diagonal(y,:,i,1,1)~=0;
    p2=Longarina_Diagonal(y,:,i,2,1)~=0;
    p=plot(Longarina_Diagonal(x,p1,i,1,1),EY*Longarina_Diagonal(y,p1,i,1,1), 'g');
    p.LineWidth = espessura_da_linha;
    p=plot(Longarina_Diagonal(x,p2,i,2,1),EY*Longarina_Diagonal(y,p2,i,2,1), 'g');
    p.LineWidth = espessura_da_linha;
end

for i=1:size(Longarina_Diagonal,3)
    p1=Longarina_Diagonal(y,:,i,1,2)~=0;
    p2=Longarina_Diagonal(y,:,i,2,2)~=0;
    p=plot(Longarina_Diagonal(x,p1,i,1,2),EY*Longarina_Diagonal(y,p1,i,1,2), '--r');
    p.LineWidth = espessura_da_linha;
    p=plot(Longarina_Diagonal(x,p2,i,2,2),EY*Longarina_Diagonal(y,p2,i,2,2), '--r');
    p.LineWidth = espessura_da_linha;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3);
scatter(EX*X_dot,EY*Y_dot,1,EI*Z_dot)
title({'Nuvem de pontos da ponte';'com dormetação e entalhes'})
xlabel('X')
ylabel('Y')
zlabel('Z')
if(dir_encoder>0)
    axis([limitInfXplot limitSupXplot margem_Inferior_Y margem_Superior_Y ])
else
    axis([limitInfXplot limitSupXplot  margem_Superior_Y margem_Inferior_Y])
end

if(EY>0)
    axis ij
else
    axis xy
end
xlabel({'X  [mm]';'C'});
ylabel('Y  [mm]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
clear m yblob;

%% contorno dos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,1,di,3);dormente_pontos(x,2,di,3);dormente_pontos(x,4,di,3);dormente_pontos(x,3,di,3);dormente_pontos(x,1,di,3)] [dormente_pontos(y,1,di,3);dormente_pontos(y,2,di,3);dormente_pontos(y,4,di,3);dormente_pontos(y,3,di,3);dormente_pontos(y,1,di,3)]];
    p=plot(han(:,1),han(:,2),':k');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end

%% realce das longarinas verticais nos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,5,di,2);dormente_pontos(x,1,di,2)] [dormente_pontos(y,5,di,2);dormente_pontos(y,1,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
    han=[[dormente_pontos(x,2,di,2);dormente_pontos(x,6,di,2)] [dormente_pontos(y,2,di,2);dormente_pontos(y,6,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,7,di,2);dormente_pontos(x,3,di,2)] [dormente_pontos(y,7,di,2);dormente_pontos(y,3,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
    han=[[dormente_pontos(x,4,di,2);dormente_pontos(x,8,di,2)] [dormente_pontos(y,4,di,2);dormente_pontos(y,8,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end
for i=1:contaux1
    han=dormente_entalhe_diagonal(y,:,i,1)~=0;
    p=plot(dormente_entalhe_diagonal(x,han,i,1),dormente_entalhe_diagonal(y,han,i,1),'--m');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
end
for i=1:contaux2
    han=dormente_entalhe_diagonal(y,:,i,2)~=0;
    p=plot(dormente_entalhe_diagonal(x,han,i,2),dormente_entalhe_diagonal(y,han,i,2),'c');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot final %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'Position', get(0, 'Screensize'));
nomeFigura="Visualização das etapas de filtragem";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end
espessura_da_linha=2.5;

figure('Name','Visualização da nuvem de pontos da ponte com dormetação e entalhes reconstruídos','NumberTitle','off');
scatter(EX*X_dot_trilho,EY*Y_dot_trilho,1,EI*Z_dot_trilho)
title({'Nuvem de pontos da ponte';'com dormetação e entalhes'})
xlabel('X')
ylabel('Y')
zlabel('Z')
if(dir_encoder>0)
    axis([limitInfXplot limitSupXplot margem_Inferior_Y margem_Superior_Y ])
else
    axis([limitInfXplot limitSupXplot  margem_Superior_Y margem_Inferior_Y])
end

if(EY>0)
    axis ij
else
    axis xy
end
xlabel({'X  [mm]';'C'});
ylabel('Y  [mm]');
hold on
clear m yblob;
%% contorno dos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,1,di,3);dormente_pontos(x,2,di,3);dormente_pontos(x,4,di,3);dormente_pontos(x,3,di,3);dormente_pontos(x,1,di,3)] [dormente_pontos(y,1,di,3);dormente_pontos(y,2,di,3);dormente_pontos(y,4,di,3);dormente_pontos(y,3,di,3);dormente_pontos(y,1,di,3)]];
    p=plot(han(:,1),han(:,2),':k');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end

%% realce das longarinas verticais nos dormentes
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,5,di,2);dormente_pontos(x,1,di,2)] [dormente_pontos(y,5,di,2);dormente_pontos(y,1,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
    han=[[dormente_pontos(x,2,di,2);dormente_pontos(x,6,di,2)] [dormente_pontos(y,2,di,2);dormente_pontos(y,6,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end
for di=1:size(dormente_salvo,3)
    han=[[dormente_pontos(x,7,di,2);dormente_pontos(x,3,di,2)] [dormente_pontos(y,7,di,2);dormente_pontos(y,3,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
    han=[[dormente_pontos(x,4,di,2);dormente_pontos(x,8,di,2)] [dormente_pontos(y,4,di,2);dormente_pontos(y,8,di,2)]];
    p=plot(han(:,1),han(:,2),'-.b');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
    clear han
end
for i=1:contaux1
    han=dormente_entalhe_diagonal(y,:,i,1)~=0;
    p=plot(dormente_entalhe_diagonal(x,han,i,1),dormente_entalhe_diagonal(y,han,i,1),'--m');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
end
for i=1:contaux2
    han=dormente_entalhe_diagonal(y,:,i,2)~=0;
    p=plot(dormente_entalhe_diagonal(x,han,i,2),dormente_entalhe_diagonal(y,han,i,2),'c');
    if(size(p,1)>0)
        p.LineWidth = espessura_da_linha;
    end
end
set(gcf, 'Position', get(0, 'Screensize'));

nomeFigura="Visualização da nuvem de pontos da ponte com dormetação e entalhes reconstruídos";
if(salvarPNG)
    saveas(gcf,destinoImagens+nomeFigura+'.png')
end
if(salvarFIG)
    saveas(gcf,destinoImagens+nomeFigura+'.fig')
end
if(salvarEPS)
    saveas(gcf,destinoImagens+nomeFigura,'epsc')
end

%% salva os pontos importantes em excel para exportar pro inventor
dormente_pontos(y,:,:,:)=-dormente_pontos(y,:,:,:);
pontos=[0,0];
j=3;
for i=1:size(dormente_pontos,3)
    han=dormente_pontos(y,:,i,j)~=0;
    pontos=[pontos;[dormente_pontos(x,han,i,j);EY*dormente_pontos(y,han,i,j)]'];
end

delete dormente_dimensões.xls
pontos(:,x)=pontos(:,1)*10;%*1.004982; %% ajuste de dimensão
pontos(:,y)=pontos(:,2)*10;%*0.983792;
xlswrite('dormente_dimensões', pontos(2:end,:));

pontos=[0,0];
j=2;
for i=1:size(dormente_pontos,3)
    han=dormente_pontos(y,:,i,j)~=0;
    pontos=[pontos;[dormente_pontos(x,han,i,j);EY*dormente_pontos(y,han,i,j)]'];
end

delete dormente_entalhe_vertical.xls
pontos(:,x)=pontos(:,1)*10; %% ajuste de dimensão
pontos(:,y)=pontos(:,2)*10;
xlswrite('dormente_entalhe_vertical', pontos(2:end,:));

pontos=[0,0];
for j=1:3:4
    for i=1:size(dormente_pontos,3)
        han=dormente_pontos(y,:,i,j)~=0;
        pontos=[pontos;[dormente_pontos(x,han,i,j);EY*dormente_pontos(y,han,i,j)]'];
    end
end

delete dormente_entalhe_diagonal_e_contraventamentos.xls
pontos(:,x)=pontos(:,1)*10;%*1.004982; %% ajuste de dimensão
pontos(:,y)=pontos(:,2)*10;%*0.983792;
xlswrite('dormente_entalhe_diagonal_e_contraventamentos', pontos(2:end,:));

%xlswrite('NomeDoArquivoXLS', todos_pontos_dormentes);
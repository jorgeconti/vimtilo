tic;
disp("Carregando arquivos");
laser1Encoder=laser1Encoder(1:reducao:end,1:reducao:1280);disp("8 laser1Encoder "+toc+ "Seg");tic
laser1X=laser1X(1:reducao:end,1:reducao:1280);disp("7 laser1X "+toc+ "Seg");tic;
laser1Z=laser1Z(1:reducao:end,1:reducao:1280);disp("6 laser1Z "+toc+ "Seg");tic;
laser1Intensidade=laser1Intensidade(1:reducao:end,1:reducao:1280);disp("5.5 laser1Intensidade "+toc+ "Seg");tic;
laser1Amostragem=laser1Amostragem(1:reducao:end,:);disp("5 laser1Amostragem "+toc+ "Seg");tic;
laser2Encoder=laser2Encoder(1:reducao:end,1:reducao:1280);disp("4 laser2Encoder "+toc+ "Seg");tic;
laser2X=laser2X(1:reducao:end,1:reducao:1280);disp("3 laser2X "+toc+ "Seg");tic;
laser2Z=laser2Z(1:reducao:end,1:reducao:1280);disp("2 laser2Z "+toc+ "Seg");tic;
laser2Intensidade=laser2Intensidade(1:reducao:end,1:reducao:1280);disp("1.5 laser2Intensidade "+toc+ "Seg");tic;
laser2Amostragem=laser2Amostragem(1:reducao:end,:);disp("1 laser2Amostragem "+toc+ "Seg");tic;
disp("Arquivos carregados");


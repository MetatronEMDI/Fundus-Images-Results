%% 1. Elección y muestra de la imagen a analizar
clear
close all
clc
 
warning off

[filename, pathname] = uigetfile({'*.jpg'; '*.png'; '*.tif'; '*.ppm'}, 'File Selector');
fullpathname = strcat(pathname, filename);
Test_Image = imread(fullpathname);
figure(1); imshow(Test_Image);

R = Test_Image;
G = Test_Image;
B = Test_Image;

R(:,:,2:3)=0;
G(:,:,1)=0;
G(:,:,3)=0;
B(:,:,1:2)=0;

% figure();imshow(R);
% figure();imshow(G);
% figure();imshow(B);
%% 2. Pre-procesamiento de la imagen a analizar
Resized_Image = imresize(Test_Image, [600 600]);
% Resized_Image = imresize(R, [600 600]);
% Resized_Image = imresize(G, [600 600]);
% Resized_Image = imresize(B, [600 600]);
Converted_Image = im2double(Resized_Image); 
Lab_Image = rgb2lab(Converted_Image); 
fill = cat(3, 1,0,0); 
Filled_Image = bsxfun(@times, fill, Lab_Image); 
Reshaped_Lab_Image = reshape(Filled_Image, [], 3); 
[C, S] = pca(Reshaped_Lab_Image); 
S = reshape(S, size(Lab_Image));
S = S(:, :, 1);
Gray_Image = (S-min(S(:)))./(max(S(:))-min(S(:))); 
% figure();imshow(Gray_Image)
Enhanced_Image = adapthisteq(Gray_Image, 'numTiles', [8 8], 'nBins', 128);
% Enhanced_Image = adapthisteq(Gray_Image, 'clipLimit',0.02);
Avg_Filter = fspecial('average', [9 9]); 
Filtered_Image = imfilter(Enhanced_Image, Avg_Filter); 
Substracted_Image = imsubtract(Filtered_Image, Enhanced_Image);
% Substracted_Image = imsubtract(Enhanced_Image, Filtered_Image);
% figure();imshow(Substracted_Image);
%% 3. Umbralización de la imagen
level = Threshold_Level(Substracted_Image);
Binary_Image = imbinarize(Substracted_Image, level-0.008);
Clean_Image = bwareaopen(Binary_Image, 100); 
% figure(); imshow(Clean_Image)
Complemented_Image = imcomplement(Clean_Image); 

%% 4. Procesamiento de la máscara binaria
Final_Result = Colorize_Image(Resized_Image, Complemented_Image, [0 0 0]);

Gray_Image = rgb2gray(Final_Result);

s=strel('disk',3); 
F=imsubtract(imadd(Gray_Image,imtophat(Gray_Image,s)),imbothat(Gray_Image,s));

Adjust=imadjust(F); 

sigma=0.3;
alpha=0.1;
 
num_niveles=16;
 
Laplacian_Filter = locallapfilt(Adjust, sigma, alpha, 'NumIntensityLevels', num_niveles); 
Laplacian_Filter1 = Laplacian_Filter(:);

%% 5. Creación de mapas binaries y cálculo de estimación de saturación de oxígeno
Order = sort(Laplacian_Filter1);

Order = Order';
 
mapa_completo = Order >= 1;
num_filtrados = Order(mapa_completo);
num_filtrados_filas = num_filtrados';

% Cálculo de la saturación de Hemoglobina desoxigenada (DHb)
mapa_DHb = Order >= 1 & Order <= 95;
num_DHb = Order(mapa_DHb);
num_DHb_filas = num_DHb'; 
NGO_DHb = (num_DHb./102)*40;
sumatoria_DHb = sum(sum(NGO_DHb));
n_DHb = length(num_DHb);
SatO_DHb = sumatoria_DHb/n_DHb;
 
% Cálculo de la saturación de oxígeno hemoglobina (Hb)
mapa_Hb = Order >= 96 & Order <= 199;
num_Hb = Order(mapa_Hb);
num_Hb_filas = num_Hb';
NGO_Hb = (num_Hb./199)*53;
sumatoria_Hb = sum(sum(NGO_Hb));
n_Hb = length(num_Hb);
SatO_Hb = sumatoria_Hb/n_Hb;
 
% Cálculo de la saturación de oxígeno oxihemoglobina (HbO2)
mapa_HbO2 = Order >= 200 & Order <= 255;
num_HbO2 = Order(mapa_HbO2);
num_HbO2_filas = num_HbO2';
NGO_HbO2 = (num_HbO2./255)*99;
sumatoria_HbO2 = sum(sum(NGO_HbO2));
n_HbO2 = length(num_HbO2);
SatO_HbO2 = sumatoria_HbO2/n_HbO2;
 
% Cálculo de la saturación de oxígeno total

SO_HbO2 = (sumatoria_Hb + sumatoria_HbO2);
SO_HbHbO2 = (sumatoria_DHb + sumatoria_Hb + sumatoria_HbO2);
 
SatO = (SO_HbO2 / SO_HbHbO2) * 100;


%% 6. Asignación colorimétrica a lo propuesto como DHb, Hb, HbO2 y total de SatO2
[p,q,~]=size(Laplacian_Filter);
for i=1:1:p
    for j=1:1:q
        if (Laplacian_Filter(i,j)< 1)
            % Color Negro
            x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 1) && (Laplacian_Filter(i,j)< 32)
            % Color Morado
            x(i,j,1)=170;
            x(i,j,2)=0;
            x(i,j,3)=255; 
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 33) && (Laplacian_Filter(i,j)< 64)
            % Color Violeta
            x(i,j,1)=85;
            x(i,j,2)=0;
            x(i,j,3)=255;
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 65) && (Laplacian_Filter(i,j)< 95)
            % Color Azul
            x(i,j,1:2)=0;
            x(i,j,3)=255;
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 96) && (Laplacian_Filter(i,j)< 127)
            % Color Cyan
            x(i,j,1)=0;
            x(i,j,2:3)=255;
        elseif (Laplacian_Filter(i,j)>= 128) && (Laplacian_Filter(i,j)< 159)
            % Color Verde
            x(i,j,1)=0;
            x(i,j,2)=187;
            x(i,j,3)=0;
        elseif (Laplacian_Filter(i,j)>= 160) && (Laplacian_Filter(i,j)< 191)
            % Color Amarillo 
            x(i,j,1:2)=255;
            x(i,j,3)=0;
        elseif (Laplacian_Filter(i,j)>= 192) && (Laplacian_Filter(i,j)< 223)
            % Color Naranja
            x(i,j,1)=255;
            x(i,j,2)=136;
            x(i,j,3)=0;
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 224) && (Laplacian_Filter(i,j)< 255)
            % Color Rojo
            x(i,j,1)=255;
            x(i,j,2:3)=0;    
%             x(i,j,1:2:3)=0;
        end
    end
end
x=x/255;
% figure();imshow(x)

% Asignación de colores Desoxihemoglobina (DHb)
[p,q,~]=size(Laplacian_Filter);
for i=1:1:p
    for j=1:1:q
        if (Laplacian_Filter(i,j)< 1)
            % Color Negro
            x1(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 1) && (Laplacian_Filter(i,j)< 64)
            % Color Morado
            x1(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 65) && (Laplacian_Filter(i,j)< 95)
            % Color Azul
            x1(i,j,1:2)=0;
            x1(i,j,3)=255;
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 96) && (Laplacian_Filter(i,j)< 255)
            x1(i,j,1:2:3)=0;
        end
    end
end
x1=x1/255;
 
% Asignación de colores Oxihemoglobina venosa (HbO2_ven)
[p,q,~]=size(Laplacian_Filter);
for i=1:1:p
    for j=1:1:q
        if (Laplacian_Filter(i,j)< 1)
            % Color Negro
            x2(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 1) && (Laplacian_Filter(i,j)< 95)
            % Color Morado
            x2(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 96) && (Laplacian_Filter(i,j)< 127)
            % Color Cyan
            x2(i,j,1)=0;
            x2(i,j,2:3)=255;
        elseif (Laplacian_Filter(i,j)>= 128) && (Laplacian_Filter(i,j)< 159)
            % Color Verde
            x2(i,j,1)=0;
            x2(i,j,2)=187;
            x2(i,j,3)=0;
        elseif (Laplacian_Filter(i,j)>= 160) && (Laplacian_Filter(i,j)< 191)
            % Color Amarillo 
            x2(i,j,1:2)=255;
            x2(i,j,3)=0;
        elseif (Laplacian_Filter(i,j)>= 192) && (Laplacian_Filter(i,j)< 255)
            x2(i,j,1:2:3)=0;
        end
    end
end
x2=x2/255;
 
% Asignación de colores Oxihemoglobina (HbO2)
[p,q,r]=size(Laplacian_Filter);
for i=1:1:p
    for j=1:1:q
        if (Laplacian_Filter(i,j)< 1)
            % Color Negro
            x3(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 1) && (Laplacian_Filter(i,j)< 191)
            % Color Negro
            x3(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 192) && (Laplacian_Filter(i,j)< 223)
            % Color Naranja
            x3(i,j,1)=255;
            x3(i,j,2)=136;
            x3(i,j,3)=0;
%             x(i,j,1:2:3)=0;
        elseif (Laplacian_Filter(i,j)>= 224) && (Laplacian_Filter(i,j)< 255)
            % Color Rojo
            x3(i,j,1)=255;
            x3(i,j,2:3)=0;        
%             x(i,j,1:2:3)=0;
        end
    end
end
x3=x3/255;


%% 7. Barra de color y muestra de las imágenes resultantes
CMap = [170,  0, 255         %// purple
        85,   0, 255         %// violet
         0,   0, 255         %// Blue
         0, 255, 255         %// Cyan
         0, 187,   0         %// Green
        255, 255,  0         %// Yellow
        255, 136,  0         %// Orange
        255,  0,   0]./255;  %// Red
        
        
Interpolate = [0
               36.42
               72.85
               109.28
               145.71
               182.44
               219.17
               255];
 
map = interp1(Interpolate/255,CMap,linspace(0,1,255));

% Saturación de oxígeno total
figure(6), imshow(x);
caxis([0 100]);
c=colorbar;
c.Label.String=['SatO2: ' num2str(round(SatO)) ' %'];
c.Label.FontSize = 20;
New_Colormap=map;
colormap(New_Colormap);

% Saturación de oxígeno total
figure(7), imshowpair(Resized_Image, x, 'montage')
%title('Comparación de Imágenes')
caxis([0 100]);
c=colorbar;
c.Label.String=['SatO2: ' num2str(round(SatO)) ' %'];
c.Label.FontSize = 20;
New_Colormap=map;
colormap(New_Colormap);

% Saturación de oxígeno Desoxihemoglobina (DHb)
figure(8), imshow(x1);
caxis([0 100]);
c=colorbar;
c.Label.String=['SatO2 - Hb: ' num2str(round(SatO_DHb)) ' %'];
c.Label.FontSize = 20;
New_Colormap=map;
colormap(New_Colormap);

% Saturación de oxígeno Hemoglobina (HbO2_ven)
figure(9), imshow(x2);
caxis([0 100]);
c=colorbar;
c.Label.String=['SatO2 - Hb: ' num2str(round(SatO_Hb)) ' %'];
c.Label.FontSize = 20;
New_Colormap=map;
colormap(New_Colormap);

% Saturación de oxígeno Oxihemoglobina (HbO2)
figure(10), imshow(x3);
caxis([0 100]);
c=colorbar;
c.Label.String=['SatO2 - HbO2: ' num2str(round(SatO_HbO2)) ' %'];
c.Label.FontSize = 20;
New_Colormap=map;
colormap(New_Colormap);

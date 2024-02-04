clear all;
close all;

% ================================================================
% Partie II : communications avec les mobiles par satellites     =
%                                                                =
%                                                                =
% 2 - Etude de la diversité apportée par le codage               =
% dans un canal de Rayleigh non sélectif en fréquence            =
%                                                                =
% Question 4 :                                                   =
% Cas de la modulation BPSK                                      = 
% Vérifiez qu’avec un entrelaceur adéquat on obtient             =
% les valeurs de diversité attendues pour un égaliseur ML.       =
%                                                                = 
% ================================================================

% Constantes
N = 20000;
Fe = 24000;
Rb = 3000;
Tb = 1/Rb;
SPAN = 4;

% Modulation QPSK - 4 symboles
n = 2;    % nb de bits par symboles
M = 2^n;  % span 4 symbols car QPSK

Rs = Rb/n;
Ns = Fe/Rs; % nbe échantillons / symboles
Ts = 1/Rs;
Te = 1/Fe;

% Génération aléatoire des N bits = information
bits = randi([0 1], 1, N);

% Entrelaceur
bits_code = convintrlv(bits,8,3);

% Code convoutif (3,1/2)
% Définition des polynomes générateurs
g1 = 5; % 5 octal -> 101
g2 = 7; % 7 octal -> 111

% Création du treillis pour ce code convolutif (3,1/2)
trellis = poly2trellis(3, [g1 g2]); 
bits_code = convenc(bits_code, trellis);  


% Mapping des bits en symboles QPSK
% Génération des symboles
% On passe les bits en +1 -1
symboles_bits = bits_code * 2 - 1;

%  1 + j = 00 = -1 -1
% -1 + j = 10 =  1 -1
% -1 - j = 11 =  1  1
%  1 - j = 01 = -1  1

% Génération des symboles
symboles = zeros(1, length(bits)/2);
for i = 1:2:length(bits)
    symboles((i+1)/2) = -symboles_bits(i) - 1i * symboles_bits(1+i);
end


% Emetteur
signal = kron(symboles, [1 zeros(1, Ns -1)]);   % On prolonge de symbole avec N-1 zéros derrière

% Filtre de mise en forme 
% Racine de cosinus surélevé Beta = 0.35
Beta = 0.35;
h_n = rcosdesign(Beta, SPAN, Ns,"sqrt"); 

% Signal en sortie de l'émetteur
x = conv(signal, h_n, 'full');

%% Implantion du canal de Rayleigh %%%%%%%%%%
% Paramètes du canal
Tc = 1; % Temps de corrélation du canal
K = 3;  

% Canal de Rayleigh
c = K*(randn(1, length(symboles)/Tc) + 1i*randn(1, length(symboles)/Tc));
c = kron(c, [1 zeros(1, Tc -1)]); % On récupère Tc échantillons identiques

z_r = symboles .* c; % Signal après passage dans le canal de Rayleigh

% % Tracé du signal modulé
% figure;
% plot(real(x), 'b');
% hold on
% plot(imag(x), 'r');
% title('Signal modulé (partie réelle et imaginaire)');
% xlabel('Temps');
% ylabel('Amplitude');
% legend('Partie réelle', 'Partie imaginaire');
% grid on;
% hold off;

% ajout d'un canal AWGN à la chaine de transmission
Eb_sur_N0 = -8:4;   % Eb/N0 variant de 0dB à 4dB pas pas de 1dB

% Je ne pense pas que la boucle for soit encore correcte car quand on fait
% varier le SNR de 0 à 4 on ne il n'y a pas de décroissance du TEB et TES

% Création de la liste TEB pour les comparer
mat_TEB_hard = zeros(1,length(Eb_sur_N0));
mat_TEB_soft = zeros(1,length(Eb_sur_N0));


for k = 1:length(Eb_sur_N0)

    % Calcul de la puissance du bruit en fonction du SNR
    Eb_sur_N0_lineaire= 10^(Eb_sur_N0(k) / 10); % Conversion du SNR en échelle linéaire
    
    % Calculez la puissance du signal à partir du signal modulé (x)
    variance_bruit = (var(x)*Ns)/(2*log2(M)*Eb_sur_N0_lineaire);
    
    % Générez un bruit gaussien pour les parties réelles et imaginaires
    bruit_reel = sqrt(variance_bruit / 2) * randn(size(z_r));
    bruit_imag = sqrt(variance_bruit / 2) * randn(size(z_r));



    % signal à l'entrée du récepteur
    x_bruite = z_r + bruit_reel + 1i*bruit_imag;


    % Récepteur

    % Egalisateur ML
    h_ml = conj(c);
    x_ml = x_bruite .* h_ml;

    
    % Demapping

    % ML
    % Décision sur les parties réelles et imaginaires
    % Reconstruction des symboles estimés à partir des décisions R et I
    symboles_estimes_ml = ((real(x_ml) > 0)*2-1) + 1i *((imag(x_ml) > 0)*2-1); 
   

    % Décision sur la partie imaginaire et réel direct
    symbole_recu_ml = zeros(1,2*length(symboles_estimes_ml));
    for i = 1:2:2*length(symboles_estimes_ml)
        symbole_recu_ml(i) = -real(symboles_estimes_ml((i+1)/2));
        symbole_recu_ml(i+1) = -imag(symboles_estimes_ml((i+1)/2));
    end
    
    bits_recu = (symbole_recu_ml + 1) /2;


    % Hard
    % On se peut donc maitenant reconstituer l'information envoyé
    % Décodage du code convolutif avec l'algorithme de Viterbi Hard
    bits_decode_hard = vitdec(bits_recu, trellis, 5, 'trunc', 'hard');


    % Soft
    % Décodage du code convolutif avec l'algorithme de Viterbi Soft
    % On crée notre vecteur alterné reel - imag
    Vect_soft = zeros(1, 2*length(x_ml));
    Vect_soft(1:2:end) = real(x_ml);
    Vect_soft(2:2:end) = imag(x_ml);

    % Décodage Viterbi Soft
    bits_decode_soft = vitdec(Vect_soft, trellis, 5, 'trunc', 'unquant');


    % Desentrelaceur
    % HARD
    bits_hard = convdeintrlv(bits_decode_hard, 8, 3);

    % SOFT
    bits_soft = convdeintrlv(bits_decode_soft, 8, 3);


    % Calcul TEB 
    % Ajoute de TEB dans la matrice
    mat_TEB_ml_soft(k) = sum(bits_soft ~= bits) / N;
    mat_TEB_ml_hard(k) = sum(bits_hard ~= bits) / N;
    
end

% Tracer TES et TEB sur le meme graphique
figure;
plot(Eb_sur_N0, mat_TEB_ml_soft, 'b', 'DisplayName', 'TEB_ML_soft');
hold on;
plot(Eb_sur_N0, mat_TEB_ml_hard, 'r', 'DisplayName', 'TEB_ML_hard');
hold off;

title("Taux d'erreur en fonction du SNR");
xlabel('SNR (dB)');
ylabel('Taux d''erreur');
legend('TEB ML soft','TEB ML hard');
grid on;

% Calcul de la diversité
p_ml_soft = polyfit(10.^(Eb_sur_N0/10), log(mat_TEB_ml_soft), 1);
p_ml_hard = polyfit(10.^(Eb_sur_N0/10), log(mat_TEB_ml_hard), 1);

diversite_ml_soft = -p_ml_soft(1);
diversite_ml_hard = -p_ml_hard(1);

% Affichage de la diversité
disp("Diversité ML soft : " + diversite_ml_soft);
disp("Diversité ML hard : " + diversite_ml_hard);
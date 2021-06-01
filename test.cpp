//
//  main.cpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//
#include <BSplineCPPConfig.h>

#include <iostream>
#include "BSpline.hpp"
#include "Parametizer.hpp"

float cps[] =
{
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.02787076249348602, -0.017719376641930098, -0.010315414918655091, 0.014146403438079295, -0.006889147881843044, -0.012896169263600172, 0.006283828771088784,
  0.05574152498697204, -0.035438753283860196, -0.020630829837310183, 0.02829280687615859, -0.013778295763686088, -0.025792338527200344, 0.012567657542177568,
  0.08361228748045806, -0.0531581299257903, -0.030946244755965274, 0.042439210314237884, -0.02066744364552913, -0.03868850779080051, 0.018851486313266354,
  0.11148304997394408, -0.07087750656772039, -0.041261659674620366, 0.05658561375231718, -0.027556591527372176, -0.05158467705440069, 0.025135315084355136,
  0.1393538124674301, -0.08859688320965049, -0.05157707459327546, 0.07073201719039647, -0.03444573940921522, -0.06448084631800086, 0.03141914385544392,
  0.1672245749609161, -0.1063162598515806, -0.06189248951193055, 0.08487842062847577, -0.04133488729105826, -0.07737701558160102, 0.03770297262653271,
  0.19509533745440216, -0.1240356364935107, -0.07220790443058565, 0.09902482406655508, -0.04822403517290131, -0.09027318484520121, 0.0439868013976215,
  0.22296609994788816, -0.14175501313544078, -0.08252331934924073, 0.11317122750463436, -0.05511318305474435, -0.10316935410880138, 0.05027063016871027,
  0.2508368624413742, -0.15947438977737088, -0.09283873426789584, 0.12731763094271367, -0.0620023309365874, -0.11606552337240156, 0.05655445893979906,
  0.2787076249348602, -0.17719376641930099, -0.10315414918655091, 0.14146403438079294, -0.06889147881843044, -0.1289616926360017, 0.06283828771088784,
  0.30657838742834626, -0.1949131430612311, -0.11346956410520602, 0.15561043781887227, -0.07578062670027348, -0.1418578618996019, 0.06912211648197664,
  0.3344491499218322, -0.2126325197031612, -0.1237849790238611, 0.16975684125695154, -0.08266977458211652, -0.15475403116320205, 0.07540594525306542,
  0.36231991241531825, -0.23035189634509126, -0.1341003939425162, 0.18390324469503083, -0.08955892246395956, -0.16765020042680223, 0.08168977402415419,
  0.3901906749088043, -0.2480712729870214, -0.1444158088611713, 0.19804964813311016, -0.09644807034580262, -0.18054636969040241, 0.087973602795243,
  0.4180614374022903, -0.2657906496289515, -0.15473122377982637, 0.21219605157118943, -0.10333721822764566, -0.1934425389540026, 0.09425743156633176,
  0.4459321998957763, -0.28351002627088157, -0.16504663869848146, 0.22634245500926872, -0.1102263661094887, -0.20633870821760275, 0.10054126033742054,
  0.47380296238926234, -0.30122940291281164, -0.17536205361713655, 0.24048885844734802, -0.11711551399133174, -0.2192348774812029, 0.10682508910850932,
  0.5016737248827484, -0.31894877955474177, -0.18567746853579167, 0.25463526188542734, -0.1240046618731748, -0.23213104674480312, 0.11310891787959812,
  0.5295444873762344, -0.3366681561966719, -0.19599288345444676, 0.2687816653235066, -0.13089380975501783, -0.24502721600840327, 0.1193927466506869,
  0.5574152498697204, -0.35438753283860197, -0.20630829837310183, 0.2829280687615859, -0.13778295763686088, -0.2579233852720034, 0.12567657542177568,
  0.5852860123632064, -0.37210690948053204, -0.21662371329175692, 0.2970744721996652, -0.1446721055187039, -0.2708195545356036, 0.13196040419286448,
  0.6131567748566925, -0.3898262861224622, -0.22693912821041204, 0.31122087563774453, -0.15156125340054696, -0.2837157237992038, 0.13824423296395327,
  0.6410275373501785, -0.40754566276439225, -0.23725454312906713, 0.3253672790758238, -0.15845040128239002, -0.29661189306280394, 0.14452806173504204,
  0.6688982998436644, -0.4252650394063224, -0.2475699580477222, 0.33951368251390307, -0.16533954916423305, -0.3095080623264041, 0.15081189050613084,
  0.6967690623371504, -0.44298441604825245, -0.2578853729663773, 0.35366008595198234, -0.17222869704607607, -0.3224042315900043, 0.1570957192772196,
  0.7246398248306365, -0.4607037926901825, -0.2682007878850324, 0.36780648939006166, -0.17911784492791913, -0.33530040085360446, 0.16337954804830837,
  0.7525105873241226, -0.47842316933211265, -0.2785162028036875, 0.381952892828141, -0.18600699280976218, -0.3481965701172047, 0.16966337681939717,
  0.7803813498176086, -0.4961425459740428, -0.2888316177223426, 0.3960992962662203, -0.19289614069160524, -0.36109273938080483, 0.175947205590486,
  0.8082521123110945, -0.5138619226159729, -0.2991470326409977, 0.4102456997042995, -0.19978528857344827, -0.373988908644405, 0.18223103436157473,
  0.8361228748045806, -0.531581299257903, -0.30946244755965274, 0.42439210314237885, -0.20667443645529132, -0.3868850779080052, 0.18851486313266352,
  0.8639936372980666, -0.549300675899833, -0.3197778624783078, 0.4385385065804581, -0.21356358433713435, -0.3997812471716053, 0.1947986919037523,
  0.8918643997915526, -0.5670200525417631, -0.3300932773969629, 0.45268491001853745, -0.2204527322189774, -0.4126774164352055, 0.2010825206748411,
  0.9197351622850387, -0.5847394291836933, -0.34040869231561804, 0.46683131345661677, -0.22734188010082046, -0.4255735856988057, 0.20736634944592988,
  0.9476059247785247, -0.6024588058256233, -0.3507241072342731, 0.48097771689469604, -0.2342310279826635, -0.4384697549624058, 0.21365017821701865,
  0.9754766872720108, -0.6201781824675534, -0.3610395221529282, 0.49512412033277536, -0.24112017586450654, -0.451365924226006, 0.21993400698810744,
  1.0033474497654968, -0.6378975591094835, -0.37135493707158335, 0.5092705237708547, -0.2480093237463496, -0.46426209348960623, 0.22621783575919624,
  1.0312182122589828, -0.6556169357514137, -0.3816703519902384, 0.523416927208934, -0.2548984716281926, -0.47715826275320633, 0.232501664530285,
  1.0590889747524688, -0.6733363123933438, -0.39198576690889353, 0.5375633306470132, -0.26178761951003565, -0.49005443201680654, 0.2387854933013738,
  1.0869597372459547, -0.6910556890352738, -0.40230118182754854, 0.5517097340850925, -0.2686767673918787, -0.5029506012804067, 0.24506932207246257,
  1.1148304997394407, -0.7087750656772039, -0.41261659674620366, 0.5658561375231718, -0.27556591527372176, -0.5158467705440068, 0.25135315084355137,
  1.142701262232927, -0.7264944423191341, -0.4229320116648588, 0.5800025409612511, -0.2824550631555648, -0.5287429398076071, 0.25763697961464016,
  1.1705720247264129, -0.7442138189610641, -0.43324742658351384, 0.5941489443993304, -0.2893442110374078, -0.5416391090712072, 0.26392080838572896,
  1.1984427872198988, -0.7619331956029942, -0.44356284150216896, 0.6082953478374097, -0.2962333589192509, -0.5545352783348074, 0.27020463715681775,
  1.226313549713385, -0.7796525722449243, -0.4538782564208241, 0.6224417512754891, -0.3031225068010939, -0.5674314475984076, 0.27648846592790655,
  1.2541843122068708, -0.7973719488868544, -0.46419367133947914, 0.6365881547135682, -0.31001165468293695, -0.5803276168620077, 0.2827722946989953,
  1.282055074700357, -0.8150913255287845, -0.47450908625813426, 0.6507345581516476, -0.31690080256478004, -0.5932237861256079, 0.2890561234700841,
  1.309925837193843, -0.8328107021707146, -0.48482450117678927, 0.6648809615897269, -0.323789950446623, -0.606119955389208, 0.2953399522411728,
  1.337796599687329, -0.8505300788126448, -0.4951399160954444, 0.6790273650278061, -0.3306790983284661, -0.6190161246528082, 0.30162378101226167,
  1.365667362180815, -0.8682494554545749, -0.5054553310140996, 0.6931737684658855, -0.3375682462103092, -0.6319122939164085, 0.30790760978335047,
  1.3935381246743008, -0.8859688320965049, -0.5157707459327546, 0.7073201719039647, -0.34445739409215215, -0.6448084631800086, 0.3141914385544392,
  1.421408887167787, -0.903688208738435, -0.5260861608514097, 0.7214665753420441, -0.35134654197399523, -0.6577046324436088, 0.320475267325528
};

float knots[100];

int main(int argc, const char * argv[])
{
    int cpCount = sizeof(cps) / sizeof(float) / 7;
    BSpline spline(cps, knots, cpCount);
    spline.init(7, cpCount);

    Parametizer param(spline);
    param.init();

    vector<float> parametization = param.parametizeSigmoidal(10);

    cout << "Spline length: " << param.length << endl;
    for(int i = 1; i < parametization.size(); i++) {
        float t0 = parametization[i-1];
        float t1 = parametization[i];
        float a0 = param.arcLength(t0);
        float a1 = param.arcLength(t1);
        cout << "Arc " << i << " t: " << t1 << " l: " << a1 << " d: " << a1 - a0 << endl;
    }

    float goal[7];
    for(int j = 0; j < 10; j++) {
        spline.eval(parametization[j], goal);
        cout << " (" << goal[0];
        for(int i = 1; i < 7; i++) {
            cout << ", " << goal[i];
        }
        cout << ")" << endl;
    }

    return 0;
}

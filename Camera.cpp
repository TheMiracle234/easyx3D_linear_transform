#include "Camera.h"

#include "Surface.h"
#include "Pixel.h"

static bool change_q = false;
static bool change_p = false;
static bool change_e = false;

void Camera::showMuis()
{
    double now_mui[3][3] = { {1,0,0},{0,1,0},{0,0,1}, };
    spanChange(now_mui);
    Change(now_mui);

    char str[4][30];
    sprintf_s(str[0], "线性变换（相对于初始状态）：");
    sprintf_s(str[1], "/ %.2f, %.2f, %.2f \\", now_mui[0][0], now_mui[1][0], now_mui[2][0]);
    sprintf_s(str[2], "| %.2f, %.2f, %.2f |", now_mui[0][1], now_mui[1][1], now_mui[2][1]);
    sprintf_s(str[3], "\\ %.2f, %.2f, %.2f /", now_mui[0][2], now_mui[1][2], now_mui[2][2]);
    const char dis = 15;

    for (int i = 0;i < 4;++i) {
        outtextxy(0, i*dis, str[i]);
    }
}

void Camera::RenewTheta()
{

}
void Camera::RenewVec()
{

    double phi1 = theta[0][1] + THETA_B;
    double phi2 = theta[0][1] - THETA_B;
    vec[0][0] = cos(theta[0][1]) * cos(theta[0][0]);
    vec[0][1] = cos(theta[0][1]) * sin(theta[0][0]);
    vec[0][2] = sin(theta[0][1]);

    vec[1][0] = -DST * sin(theta[0][0]) + cos(theta[0][0]) * sqrt(1 - DST * DST) * cos(phi1);
    vec[1][1] = DST * cos(theta[0][0]) + sin(theta[0][0]) * sqrt(1 - DST * DST) * cos(phi1);
    vec[1][2] = +sqrt(1 - DST * DST) * sin(phi1);

    vec[2][0] = vec[1][0] + 2 * DST * sin(theta[0][0]);
    vec[2][1] = vec[1][1] - 2 * DST * cos(theta[0][0]);
    vec[2][2] = vec[1][2];

    vec[3][0] = -DST * sin(theta[0][0]) + cos(theta[0][0]) * sqrt(1 - DST * DST) * cos(phi2);
    vec[3][1] = DST * cos(theta[0][0]) + sin(theta[0][0]) * sqrt(1 - DST * DST) * cos(phi2);
    vec[3][2] = +sqrt(1 - DST * DST) * sin(phi2);

    vec[4][0] = vec[3][0] + 2 * DST * sin(theta[0][0]);
    vec[4][1] = vec[3][1] - 2 * DST * cos(theta[0][0]);
    vec[4][2] = vec[3][2];
}
void Camera::RenewSrc()
{
    for (int i = 0; i < 3; ++i)
    {
        src_RC[0][i] = (vec[2][i] - vec[1][i]) / ROW;
    }
    for (int i = 0; i < 3; ++i)
    {
        src_RC[1][i] = (vec[3][i] - vec[1][i]) / COL;
    }
}
//默认构造函数
Camera::Camera()
{
    loc[0] = 15;
    loc[1] = 5;
    loc[2] = 10;
    theta[0][0] = 3.330088;
    theta[0][1] = -0.644026;
    RenewVec();
    RenewSrc();
    pixels.clear();
}
Camera::~Camera() {
}
void Camera::spanChange(double p[3][3]){
    for (int i = 0;i < 3;++i) {
        for (int j = 0;j < 3;++j) {
            int index1 = (j + 1) % 3;
            int index2 = (j + 2) % 3;
            double tmp[3] = { p[i][0],p[i][1], p[i][2], };
            p[i][index1] = tmp[index1] * cos(alpha[j]) - tmp[index2] * sin(alpha[j]);
            p[i][index2] = tmp[index1] * sin(alpha[j]) + tmp[index2] * cos(alpha[j]);
        }
    }
}
void Camera::Change(double p[3][3])
{
    double tmp[3][3] = {
        {p[0][0],p[0][1],p[0][2],},
        {p[1][0],p[1][1],p[1][2],},
        {p[2][0],p[2][1],p[2][2],},
    };
    for (int i = 0;i < 3;++i) {
        for (int j = 0;j < 3;++j) {
            p[i][j] = tmp[i][0] * mui[0][j] + tmp[i][1] * mui[1][j] + tmp[i][2] * mui[2][j];
        }
    }
}
//vec0=-0.797812, vec1=-0.508502, vec2=-0.32391
//vec0=-0.804964, vec1=-0.509866, vec2=0.303429
void Camera::PrintInfo()
{
    std::cout << "loc:" << std::endl;
    std::cout << "\t";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << loc[i] << " ";
    }
    std::cout << "\nvec:" << std::endl;
    for (int i = 0; i < 5; ++i)
    {
        std::cout << "\t" << i << ": ";
        for (int j = 0; j < 3; ++j)
        {
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\ntheta:" << std::endl;
    for (int i = 0; i < 5; ++i)
    {
        std::cout << "\t" << i << ": ";
        for (int j = 0; j < 2; ++j)
        {
            std::cout << theta[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\nsrc:" << std::endl;
    std::cout << "\tsrc_R: ";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << src_RC[0][i] << " ";
    }
    std::cout << "\n\tsrc_C: ";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << src_RC[1][i] << " ";
    }
}

//射一次，对点的
bool Camera::shoot(Pixel& pix, int* the_screen)
{
    double cam[3] =
    {
        loc[0] + vec[1][0],
        loc[1] + vec[1][1],
        loc[2] + vec[1][2]
    };
    double t1 = 0, t2 = 0, t = 0;
    double the_dis = 0;
    for (int i = 0; i < 3; ++i)
    {
        t1 += (cam[i] - loc[i]) * vec[0][i];
        t2 += vec[0][i] * (pix.loc[i] - loc[i]);
        the_dis += (pix.loc[i] - loc[i]) * (pix.loc[i] - loc[i]);
    }
    t = t1 / t2;
    if (t <= 0)
    {
        return false;
    }
    the_dis = t * (the_dis);
    double Tvec[3] = {
        loc[0] + (pix.loc[0] - loc[0]) * t - cam[0],
        loc[1] + (pix.loc[1] - loc[1]) * t - cam[1],
        loc[2] + (pix.loc[2] - loc[2]) * t - cam[2],
    };
    int mui[2] =
    {
        (int)(((Tvec[0] * src_RC[1][1]) - (Tvec[1] * src_RC[1][0])) / ((src_RC[0][0] * src_RC[1][1]) - (src_RC[0][1] * src_RC[1][0]))),
        (int)(((Tvec[0] * src_RC[0][1]) - (Tvec[1] * src_RC[0][0])) / ((src_RC[1][0] * src_RC[0][1]) - (src_RC[1][1] * src_RC[0][0])))
    };
    if (the_screen != nullptr) {
        the_screen[0] = mui[0];
        the_screen[1] = mui[1];
    }
    return true;
}

//优化射面
bool Camera::shoot(Surface& surf, int) {
    double tmp_p[3][3] = {
        {surf.p[0][0],surf.p[0][1],surf.p[0][2],},
        {surf.p[1][0],surf.p[1][1],surf.p[1][2],},
        {surf.p[2][0],surf.p[2][1],surf.p[2][2],},
    };

    if (!surf.fix) {
        spanChange(tmp_p);
        Change(tmp_p);
    }

    Pixel pix[3] = { Pixel(tmp_p[0]), Pixel(tmp_p[1]), Pixel(tmp_p[2]) };
    int* the_screen[3];
    double ave_trian_loc[3] = {
        (pix[0].loc[0] + pix[1].loc[0] + pix[2].loc[0]) / 3,
        (pix[0].loc[1] + pix[1].loc[1] + pix[2].loc[1]) / 3,
        (pix[0].loc[2] + pix[1].loc[2] + pix[2].loc[2]) / 3,
    };

    double distance = 0;
    for (int i = 0;i < 3;++i) {
        double up = 0, down = 0;
        for (int j = 0;j < 3;++j) {
            up += (pix[i].loc[j] - loc[j]) * vec[0][j];
            down += vec[0][j] * vec[0][j];
        }
        up = up >= 0 ? up : -up;
        down = sqrt(down);
        double tmp_dis = up / down;
        if (i == 0) {
            distance = tmp_dis;
        }
        else if (tmp_dis < distance) {
            distance = tmp_dis;
        }
    }
    //std::cout << "distance = " << distance << std::endl;

    for (int i = 0;i < 3;++i) {
        the_screen[i] = new int[2];
        if (!shoot(pix[i], the_screen[i])) {
            for (int j = 0;j <= i;++j) {
                delete the_screen[j];
            }
            a1++;
            return false;
        }
    }

    std::vector<std::vector<int>> tmp{
        {the_screen[0][0],the_screen[0][1],surf.line},
        {the_screen[1][0],the_screen[1][1]},
        {the_screen[2][0],the_screen[2][1]}
    };
    pixels.push_back({ tmp,-distance,surf.color(ave_trian_loc[0], ave_trian_loc[1], ave_trian_loc[2]) });
    for (int i = 0;i < 3;++i) {
        delete the_screen[i];
    }
    return true;
}
//优化后的
void Camera::Show(int) {
    
    //std::cout << "a1 = " << a1 << std::endl;
    a1 = 0;

    std::sort(pixels.begin(), pixels.end(),
        [](const auto& a, const auto& b) {
            return a.distance < b.distance;
        });
    for (const auto& itr : pixels) {
        POINT tmp_loc[3]{
            {itr.loc[0][0],itr.loc[0][1]},
            {itr.loc[1][0],itr.loc[1][1]},
            {itr.loc[2][0],itr.loc[2][1]},
        };
        if (!itr.loc[0][2]) {
            setlinecolor(0xFFFFFF);
            setfillcolor(itr.color);
            if (showlines)
                fillpolygon(tmp_loc, 3);
            else
                solidpolygon(tmp_loc, 3);
        }
        else {
            setlinecolor(itr.color);
            polygon(tmp_loc, 3);
            //solidpolygon(tmp_loc, 3);
        }
    }
    pixels.clear();

    if (showmuis)
        showMuis();
}

inline void onePress(char& ch, char key,bool& flag) {
    if (ch != key)
        flag = false;
    if (flag) {
        ch = '0';
    }
    if (ch == key && !flag) {
        flag = true;
    }
}
void Camera::Command(char ch)
{
    double tmp = ((vec[0][1] * vec[0][1]) + (vec[0][0] * vec[0][0]));
    double delta_x = (delta_dst * vec[0][0]) / sqrt(tmp);
    double delta_y = (delta_dst * vec[0][1]) / sqrt(tmp);
    
    onePress(ch, 'p', change_p);
    onePress(ch, 'q', change_q);
    onePress(ch, 'e', change_e);

    bool nochosen = false;
    switch (ch)
    {
    case 'k':
    {
        theta[0][0] = fmod((theta[0][0] + delta_theta + 2 * PI), (2 * PI));
        break;
    }
    case 'h':
    {
        double tmp1 = theta[0][1] + delta_theta;
        theta[0][1] = (tmp1 >= -PI / 2 && tmp1 <= PI / 2) ? tmp1 : theta[0][1];
        break;
    }
    case 'j':
    {
        double tmp2 = theta[0][1] - delta_theta;
        theta[0][1] = (tmp2 >= -PI / 2 && tmp2 <= PI / 2) ? tmp2 : theta[0][1];
        break;
    }
    case 'l':
    {
        theta[0][0] = fmod((theta[0][0] - delta_theta + 2 * PI), (2 * PI));
        break;
    }
    case 'w'://前
        loc[0] += delta_x;
        loc[1] += delta_y;
        break;
    case 's'://后
        loc[0] -= delta_x;
        loc[1] -= delta_y;
        break;
    case 'd'://左
        loc[0] += delta_y;
        loc[1] -= delta_x;
        break;
    case 'a'://右
        loc[0] -= delta_y;
        loc[1] += delta_x;
        break;
    case '!'://上
        loc[2] += delta_dst;
        break;
    case '@'://下
        loc[2] -= delta_dst;
        break;
        //旋转变换
    case 'z':
        alpha[0] += d_theta;
        break;
    case 'x':
        alpha[1] += d_theta;
        break;
    case 'c':
        alpha[2] += d_theta;
        break;
    case 'b':
        alpha[0] += -d_theta;
        break;
    case 'n':
        alpha[1] += -d_theta;
        break;
    case 'm':
        alpha[2] += -d_theta;
        break;
        //一般线性变换
    case '1':
        mui[0][0] += d_mui * positive;
        break;
    case '2':
        mui[0][1] += d_mui * positive;
        break;
    case '3':
        mui[0][2] += d_mui * positive;
        break;
    case '4':
        mui[1][0] += d_mui * positive;
        break;
    case '5':
        mui[1][1] += d_mui * positive;
        break;
    case '6':
        mui[1][2] += d_mui * positive;
        break;
    case '7':
        mui[2][0] += d_mui * positive;
        break;
    case '8':
        mui[2][1] += d_mui * positive;
        break;
    case '9':
        mui[2][2] += d_mui * positive;
        break;
    case 'p':
        positive = (positive == 1 ? -1 : 1);
        break;
        //展示线性变换矩阵
    case 'q':
        showmuis = (showmuis ? false : true);
        break;
    case 'e':
        showlines = (showlines ? false : true);
        break;
    case 'r':
        alpha[0] = alpha[1] = alpha[2] = 0;
        mui[0][0] = mui[1][1] = mui[2][2] = 1;
        mui[0][1] = mui[0][2] = mui[1][0] = mui[1][2] = mui[2][0] = mui[2][1] = 0;
        break;
    default:
        nochosen = true;
        break;
    }
    if (!nochosen){
        printf("%c\n", ch);
    }
    RenewVec();
    RenewSrc();
}

void Camera::pushFunc(std::vector<Surface>& surf, std::function<double(double, double)> func,double the_max[3], double the_min[3], double step){
    for (double i = the_min[0];i+step <= the_max[0];i += step) {
        for (double j = the_min[1];j+step <= the_max[1];j += step) {
            double val[4] = {
                func(i,j),func(i + step,j),func(i,j + step),func(i + step,j + step)
            };
            double point1[3][3] = {
                {i,j,val[0]},
                {i + step,j,val[1]},
                {i,j + step,val[2]},
            };
            double point2[3][3] = {
                {i+step,j+step,val[3]},
                {i + step,j,val[1]},
                {i,j + step,val[2]},
            };
            int colors_[3] = {
                0xFF / (the_max[0] - the_min[0]) * (point1[0][0] - the_min[0]),
                0xFF / (the_max[1] - the_min[1]) * (point1[0][1] - the_min[1]),
                0xFF / (the_max[2] - the_min[2]) * (point1[0][2] - the_min[2]), };
            if (!(val[0] > the_max[2] || val[0] < the_min[2] || std::isnan(val[0]) ||
                val[1] > the_max[2] || val[1] < the_min[2] || std::isnan(val[1]) ||
                val[2] > the_max[2] || val[2] < the_min[2] || std::isnan(val[2]) )){
                surf.emplace_back(Surface(point1, colors_, "fixed_RGB"));

                if (!(val[3] > the_max[2] || val[3] < the_min[2] || std::isnan(val[3]))) {
                    surf.emplace_back(Surface(point2, colors_, "fixed_RGB"));
                }
            }
        }
    }
}

void Camera::Aplication(std::vector<Surface>& surf)
{
    char ch = 0;
    char juststart = 1;
    int control = 1.0;

    clock_t the_time1 = 0;
    clock_t the_time2 = 0;
    //const char CPUNUM = 8;

    do
    {
        BeginBatchDraw();
        cleardevice();

        //std::vector<std::thread> trd;
        //const char THREAD_COUNT = 8;
        //for (int t = 0; t < THREAD_COUNT; ++t) {
        //    trd.emplace_back([&, t]() {
        //        for (int i = t; i < surf.size(); i += THREAD_COUNT) {
        //            shoot(surf[i], control); // 若shoot无需锁，直接执行
        //        }
        //        });
        //}
        for (int i = 0; i < (int)surf.size(); ++i)
        {
            shoot(surf[i], control);
            //trd[i].join();
        }

        if (!juststart)
        {
            the_time1 = clock();
        }
        //std::cout << (double)(the_time1 - the_time2) / CLOCKS_PER_SEC << std::endl;
        double tmp = 0.05 - (double)(the_time1 - the_time2) / CLOCKS_PER_SEC;
        if (tmp > 0)
        {
            //preciseSleep(tmp);
        }
        Show(control);
        //Show(1.1);
        the_time2 = clock();
        ch = CheckPress();
        Command(ch);
        juststart = 0;

        FlushBatchDraw();
    } while (ch != '#');
}

//精确延迟（0.0000001s）
void preciseSleep(double seconds)
{
    LARGE_INTEGER frequency;
    LARGE_INTEGER start;
    LARGE_INTEGER current;

    // 获取性能计数器的频率
    QueryPerformanceFrequency(&frequency);
    // 获取当前性能计数器的值
    QueryPerformanceCounter(&start);

    double target = start.QuadPart + seconds * frequency.QuadPart;

    do
    {
        QueryPerformanceCounter(&current);
    } while (current.QuadPart < target);
}
//检查按键连续按下
char CheckPress()
{
    if (GetAsyncKeyState(VK_SPACE) & 0x8000)
    {
        return '!';
    }
    if (GetAsyncKeyState(VK_SHIFT) & 0x8000)
    {
        return '@';
    }
    if (GetAsyncKeyState('A') & 0x8000)
    {
        return 'a';
    }
    if (GetAsyncKeyState('D') & 0x8000)
    {
        return 'd';
    }
    if (GetAsyncKeyState('W') & 0x8000)
    {
        return 'w';
    }
    if (GetAsyncKeyState('S') & 0x8000)
    {
        return 's';
    }
    if (GetAsyncKeyState('H') & 0x8000)
    {
        return 'h';
    }
    if (GetAsyncKeyState('J') & 0x8000)
    {
        return 'j';
    }
    if (GetAsyncKeyState('K') & 0x8000)
    {
        return 'k';
    }
    if (GetAsyncKeyState('L') & 0x8000)
    {
        return 'l';
    }
    if (GetAsyncKeyState(VK_ESCAPE) & 0x8000)
    {
        return '#';
    }

    if (GetAsyncKeyState('Z') & 0x8000)
    {
        return 'z';
    }if (GetAsyncKeyState('X') & 0x8000)
    {
        return 'x';
    }if (GetAsyncKeyState('C') & 0x8000)
    {
        return 'c';
    }if (GetAsyncKeyState('B') & 0x8000)
    {
        return 'b';
    }if (GetAsyncKeyState('N') & 0x8000)
    {
        return 'n';
    }if (GetAsyncKeyState('M') & 0x8000)
    {
        return 'm';
    }

    if (GetAsyncKeyState(VK_NUMPAD1) & 0x8000)
    {
        return '1';
    }if (GetAsyncKeyState(VK_NUMPAD2) & 0x8000)
    {
        return '2';
    }if (GetAsyncKeyState(VK_NUMPAD3) & 0x8000)
    {
        return '3';
    }if (GetAsyncKeyState(VK_NUMPAD4) & 0x8000)
    {
        return '4';
    }if (GetAsyncKeyState(VK_NUMPAD5) & 0x8000)
    {
        return '5';
    }if (GetAsyncKeyState(VK_NUMPAD6) & 0x8000)
    {
        return '6';
    }if (GetAsyncKeyState(VK_NUMPAD7) & 0x8000)
    {
        return '7';
    }if (GetAsyncKeyState(VK_NUMPAD8) & 0x8000)
    {
        return '8';
    }if (GetAsyncKeyState(VK_NUMPAD9) & 0x8000)
    {
        return '9';
    }
    if ((GetAsyncKeyState('R') & 0x8000)) {
        return 'r';
    }
    //单次点击
    if ((GetAsyncKeyState('P') & 0x8000))
    {
        return 'p';
    }
    if ((GetAsyncKeyState('Q') & 0x8000))
    {
        return 'q';
    }
    if ((GetAsyncKeyState('E') & 0x8000))
    {
        return 'e';
    }
    return '0';
}

void ADD(double obj1[3], double obj2[3])
{
    for (int i = 0; i < 3; ++i)
    {
        obj1[i] += obj2[i];
    }
}


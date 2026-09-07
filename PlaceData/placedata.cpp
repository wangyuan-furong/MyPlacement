#include "placedata.h"

// 提取路径的目录部分
static string dirOf(const string &path)
{
    size_t p = path.find_last_of("/\\");
    if (p == string::npos)
        return ".";
    return path.substr(0, p);
}

// 判断是否为头部行（UCLA 头 / 注释 / 空行）
static bool isHeaderOrComment(const string &line)
{
    if (line.empty())
        return true;
    if (line[0] == '#')
        return true;
    if (line.rfind("UCLA", 0) == 0)
        return true;
    return false;
}

Module *PlaceData::getModuleByName(const string &name)
{
    map<string, Module *>::iterator it = moduleMap.find(name);
    if (it == moduleMap.end())
        return NULL;
    return it->second;
}

// 入口：解析 BookShelf 全套文件
bool PlaceData::ReadBookShelf(const string &auxPath)
{
    if (!ReadAuxFile(auxPath))
        return false;
    if (!ReadSclFile(sclFile))
        return false;
    if (!ReadNodesFile(nodesFile))
        return false;
    if (!ReadPlFile(plFile))
        return false;
    if (!ReadNetsFile(netsFile))
        return false;
    ReadWtsFile(wtsFile); // wts 可选，允许失败
    PlaceTerminalsOnBoundary();
    InitChipRegion();
    return true;
}

bool PlaceData::ReadAuxFile(const string &auxPath)
{
    basePath = dirOf(auxPath);
    ifstream in(auxPath.c_str());
    if (!in)
    {
        cerr << "[错误] 无法打开 aux 文件: " << auxPath << endl;
        return false;
    }
    string line;
    while (getline(in, line))
    {
        if (line.find("RowBasedPlacement") != string::npos)
        {
            size_t p = line.find(':');
            if (p == string::npos)
                continue;
            istringstream iss(line.substr(p + 1));
            iss >> nodesFile >> netsFile >> wtsFile >> plFile >> sclFile;
            break;
        }
    }
    // 拼接为绝对/相对完整路径（相对 aux 所在目录）
    if (!nodesFile.empty()) nodesFile = basePath + "/" + nodesFile;
    if (!netsFile.empty()) netsFile = basePath + "/" + netsFile;
    if (!wtsFile.empty()) wtsFile = basePath + "/" + wtsFile;
    if (!plFile.empty()) plFile = basePath + "/" + plFile;
    if (!sclFile.empty()) sclFile = basePath + "/" + sclFile;

    cout << "Reading AUX file: " << auxPath << endl;
    cout << "  " << nodesFile << " " << netsFile << " " << wtsFile << " " << plFile << " " << sclFile << endl;
    return true;
}

bool PlaceData::ReadSclFile(const string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cerr << "[错误] 无法打开 scl 文件: " << path << endl;
        return false;
    }
    string line;
    int numRows = 0;
    bool firstRow = true;
    float minX = 1e30f, minY = 1e30f, maxX = -1e30f, maxY = -1e30f;

    while (getline(in, line))
    {
        if (isHeaderOrComment(line))
            continue;
        if (line.rfind("NumRows", 0) == 0)
        {
            size_t p = line.find(':');
            if (p != string::npos)
                numRows = atoi(line.c_str() + p + 1);
            continue;
        }
        if (line.rfind("CoreRow", 0) == 0)
        {
            SiteRow row;
            float coordinate = 0, subrowOrigin = 0, numSites = 0;
            // 读取一个 CoreRow 块，直到 "End"
            while (getline(in, line))
            {
                if (line.rfind("End", 0) == 0)
                    break;
                if (isHeaderOrComment(line))
                    continue;
                istringstream iss(line);
                string key, colon;
                while (iss >> key)
                {
                    iss >> colon; // 吞掉 ":"
                    if (key == "Coordinate") iss >> coordinate;
                    else if (key == "Height") iss >> row.height;
                    else if (key == "Sitewidth") iss >> siteWidth;
                    else if (key == "Sitespacing") iss >> row.step;
                    else if (key == "Siteorient") { int v; iss >> v; row.orientation = (v == 1) ? N : S; }
                    else if (key == "SubrowOrigin") iss >> subrowOrigin;
                    else if (key == "NumSites") iss >> numSites;
                }
            }
            row.bottom = coordinate;
            row.start = POS_2D(subrowOrigin, coordinate);
            row.end = POS_2D(subrowOrigin + numSites * siteWidth, coordinate);
            SiteRows.push_back(row);

            // 更新 coreRegion 外包范围
            minX = min(minX, (float)subrowOrigin);
            minY = min(minY, (float)coordinate);
            maxX = max(maxX, (float)(subrowOrigin + numSites * siteWidth));
            maxY = max(maxY, (float)(coordinate + row.height));
            if (firstRow)
            {
                siteHeight = row.height;
                firstRow = false;
            }
        }
    }

    if (SiteRows.empty())
    {
        cerr << "[错误] scl 文件中没有 CoreRow: " << path << endl;
        return false;
    }
    coreRegion.ll = POS_2D(minX, minY);
    coreRegion.ur = POS_2D(maxX, maxY);
    cout << "Set core region from site info: ll " << coreRegion.ll << "  ur " << coreRegion.ur << endl;
    return true;
}

bool PlaceData::ReadNodesFile(const string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cerr << "[错误] 无法打开 nodes 文件: " << path << endl;
        return false;
    }
    string line;
    int numNodes = 0;
    while (getline(in, line))
    {
        if (isHeaderOrComment(line))
            continue;
        if (line.rfind("NumNodes", 0) == 0)
        {
            size_t p = line.find(':');
            if (p != string::npos)
                numNodes = atoi(line.c_str() + p + 1);
            continue;
        }
        if (line.rfind("NumTerminals", 0) == 0)
        {
            size_t p = line.find(':');
            if (p != string::npos)
                terminalCount = atoi(line.c_str() + p + 1);
            continue;
        }
        // 单元行：名称 宽 高 [terminal]
        istringstream iss(line);
        string name, term;
        float w, h;
        iss >> name >> w >> h >> term;

        Module *m = new Module();
        m->idx = (int)moduleCount++;
        m->name = name;
        m->width = w;
        m->height = h;
        m->area = w * h;
        m->isFixed = (term == "terminal");

        if (m->isFixed)
        {
            m->isMacro = false;
            Terminals.push_back(m);
        }
        else
        {
            // 宏单元判定：高度大于标准单元行高（site 高度）
            m->isMacro = (siteHeight > 0 && h > siteHeight + 1e-6f);
            if (m->isMacro)
                MacroCount++;
            Nodes.push_back(m);
        }
        moduleMap[name] = m;
    }

    // 预分配
    Nodes.reserve(Nodes.size());
    Terminals.reserve(Terminals.size());
    cout << "NumModules: " << moduleCount << endl;
    cout << "  NumNodes(可移动): " << Nodes.size() << "  Terminals: " << Terminals.size()
         << "  Macro: " << MacroCount << endl;
    return true;
}

bool PlaceData::ReadPlFile(const string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cerr << "[错误] 无法打开 pl 文件: " << path << endl;
        return false;
    }
    string line;
    while (getline(in, line))
    {
        if (isHeaderOrComment(line))
            continue;
        istringstream iss(line);
        string name, colon, orient;
        float x, y;
        iss >> name >> x >> y >> colon >> orient;
        Module *m = getModuleByName(name);
        if (m)
        {
            m->ll = POS_2D(x, y);
            m->orientation = parseOrientation(orient);
        }
    }
    cout << "Initialize module position with file: " << path << endl;
    return true;
}

bool PlaceData::ReadNetsFile(const string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cerr << "[错误] 无法打开 nets 文件: " << path << endl;
        return false;
    }
    string line;
    int numNets = 0, numPins = 0;
    while (getline(in, line))
    {
        if (isHeaderOrComment(line))
            continue;
        if (line.rfind("NumNets", 0) == 0)
        {
            size_t p = line.find(':');
            if (p != string::npos)
                numNets = atoi(line.c_str() + p + 1);
            continue;
        }
        if (line.rfind("NumPins", 0) == 0)
        {
            size_t p = line.find(':');
            if (p != string::npos)
                numPins = atoi(line.c_str() + p + 1);
            Nets.reserve(numNets);
            Pins.reserve(numPins);
            continue;
        }
        if (line.rfind("NetDegree", 0) == 0)
        {
            istringstream iss(line);
            string key, colon, netname;
            int degree;
            iss >> key >> colon >> degree >> netname;

            Net *net = new Net();
            net->idx = (int)Nets.size();
            net->netPins.reserve(degree);

            for (int i = 0; i < degree; i++)
            {
                if (!getline(in, line))
                    break;
                istringstream pins(line);
                string mname, pname, c2;
                float xoff, yoff;
                pins >> mname >> pname >> c2 >> xoff >> yoff;

                Pin *pin = new Pin();
                pin->idx = (int)Pins.size();
                pin->module = getModuleByName(mname);
                pin->net = net;
                pin->offset = POS_2D(xoff, yoff);
                pin->direction = (pname == "O") ? 0 : (pname == "I" ? 1 : -1);
                net->netPins.push_back(pin);
                Pins.push_back(pin);
                if (pin->module)
                    pin->module->modulePins.push_back(pin);
            }
            Nets.push_back(net);
        }
    }
    netCount = (int)Nets.size();
    pinCount = (int)Pins.size();
    cout << "Nets: " << netCount << "  Pins: " << pinCount << endl;
    return true;
}

bool PlaceData::ReadWtsFile(const string &path)
{
    ifstream in(path.c_str());
    if (!in)
        return true; // wts 可选，不存在则跳过
    string line;
    while (getline(in, line))
    {
        if (isHeaderOrComment(line))
            continue;
        istringstream iss(line);
        string netname;
        double w;
        iss >> netname >> w;
        // BookShelf 里 net 名在 .nets 中，这里简化：直接忽略（本设计不依赖 wts）
        (void)netname;
        (void)w;
    }
    return true;
}

// BookShelf 的 .pl 文件常把 terminal（I/O pad）写成 0 0，导致它们全部挤在原点，
// 二次解析的终端拉力退化为 0。这里检测到这种情况时，把 terminal 沿 core 四边均匀铺开。
void PlaceData::PlaceTerminalsOnBoundary()
{
    if (Terminals.empty() || coreRegion.area() <= 0)
        return;
    float coreW = coreRegion.ur.x - coreRegion.ll.x;
    float coreH = coreRegion.ur.y - coreRegion.ll.y;

    // 统计 terminal 当前位置跨度
    float minX = 1e30f, minY = 1e30f, maxX = -1e30f, maxY = -1e30f;
    for (size_t i = 0; i < Terminals.size(); i++)
    {
        minX = min(minX, Terminals[i]->ll.x);
        minY = min(minY, Terminals[i]->ll.y);
        maxX = max(maxX, Terminals[i]->ll.x);
        maxY = max(maxY, Terminals[i]->ll.y);
    }
    // 若 terminal 已明显铺开（跨度 > core 的 1%），说明 .pl 已指定位置，不做处理
    if ((maxX - minX) > coreW * 0.01f || (maxY - minY) > coreH * 0.01f)
        return;

    // 沿 core 四边按弧长均匀分布（中心落在边界上）
    int T = (int)Terminals.size();
    float W = coreW, H = coreH;
    float perim = 2.0f * (W + H);
    for (int i = 0; i < T; i++)
    {
        float s = (i + 0.5f) * perim / T; // terminal 中心在周长上的弧长
        float cx, cy;
        if (s < W) // 底边
        {
            cx = coreRegion.ll.x + s;
            cy = coreRegion.ll.y;
        }
        else if (s < W + H) // 右边
        {
            cx = coreRegion.ur.x;
            cy = coreRegion.ll.y + (s - W);
        }
        else if (s < 2 * W + H) // 顶边
        {
            cx = coreRegion.ur.x - (s - W - H);
            cy = coreRegion.ur.y;
        }
        else // 左边
        {
            cx = coreRegion.ll.x;
            cy = coreRegion.ur.y - (s - 2 * W - H);
        }
        Module *t = Terminals[i];
        t->ll = POS_2D(cx - t->width * 0.5f, cy - t->height * 0.5f);
    }
    cout << "  [终端铺开] " << T << " 个 terminal 已沿 core 边界均匀分布" << endl;
}

void PlaceData::InitChipRegion()
{
    float minX = coreRegion.ll.x, minY = coreRegion.ll.y;
    float maxX = coreRegion.ur.x, maxY = coreRegion.ur.y;
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        minX = min(minX, Nodes[i]->ll.x);
        minY = min(minY, Nodes[i]->ll.y);
        maxX = max(maxX, Nodes[i]->ur().x);
        maxY = max(maxY, Nodes[i]->ur().y);
    }
    for (size_t i = 0; i < Terminals.size(); i++)
    {
        minX = min(minX, Terminals[i]->ll.x);
        minY = min(minY, Terminals[i]->ll.y);
        maxX = max(maxX, Terminals[i]->ur().x);
        maxY = max(maxY, Terminals[i]->ur().y);
    }
    chipRegion.ll = POS_2D(minX, minY);
    chipRegion.ur = POS_2D(maxX, maxY);
}

void PlaceData::PrintSummary()
{
    double movableArea = 0, fixedArea = 0;
    for (size_t i = 0; i < Nodes.size(); i++)
        movableArea += (double)Nodes[i]->width * Nodes[i]->height;
    for (size_t i = 0; i < Terminals.size(); i++)
        fixedArea += (double)Terminals[i]->width * Terminals[i]->height;
    double coreArea = coreRegion.area();

    int maxDegree = 0;
    for (size_t i = 0; i < Nets.size(); i++)
        maxDegree = max(maxDegree, (int)Nets[i]->netPins.size());

    cout << endl
         << "<<<< DATABASE SUMMARY >>>>" << endl;
    cout << "  Core region: ll " << coreRegion.ll << "  ur " << coreRegion.ur << endl;
    cout << "  Row Height/Number: " << siteHeight << " / " << SiteRows.size() << endl;
    cout << "  Core Area: " << coreArea << endl;
    cout << "  Cell Area: " << movableArea << " (" << movableArea / coreArea * 100 << "%)" << endl;
    cout << "  Movable Area: " << movableArea << endl;
    cout << "  Fixed Area: " << fixedArea << " (" << fixedArea / coreArea * 100 << "%)" << endl;
    cout << "  Cell #: " << Nodes.size() << endl;
    cout << "  Object #: " << moduleCount << " (fixed: " << Terminals.size() << ") (macro: " << MacroCount << ")" << endl;
    cout << "  Net #: " << netCount << endl;
    cout << "  Max net degree: " << maxDegree << endl;
    cout << "  Pin #: " << pinCount << endl;
}

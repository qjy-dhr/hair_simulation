#include "projdyn_api.h"
#include <thread>
#include <chrono>
#include "projdyn.h"
#include <nanogui/slider.h>
#include "viewer.h"
#include <nanogui/theme.h>
using namespace nanogui;

// Global variables used by the callback functions
int projdyn_num_iterations = PROJDYN_NUM_ITS_INITIAL;
int projdyn_num_hairs = PROJDYN_NUM_HAIRS_INITIAL;
int projdyn_res = PROJDYN_RES_INITIAL;
float projdyn_ball_radius = PROJDYN_BALL_RADIUS_INITIAL;
float projdyn_hair_thickness = PROJDYN_HAIR_THICKNESS_INITIAL;
float projdyn_seg_length = PROJDYN_SEG_LENGTH_INITIAL;
ProjDyn::Simulator sim;
Eigen::MatrixXf upload_pos, upload_rods, upload_rods_tan, upload_rods_norm;
std::thread projdyn_thread;
bool projdyn_active = false;
bool default_constraints = true;

//Store the execution time for averaging
std::vector<int> exec_time;
size_t num_iter;

// 将更新的顶点位置传递给查看器
bool projdyn_setmesh(Viewer* viewer) {
	return projdyn_setmesh(viewer, false);
}

// 在查看器中加载新网格时，也要更新模拟器中的网格
bool projdyn_setmesh(Viewer* viewer, bool add_tets) {
	std::cout << "New mesh was loaded, re-initializing simulation..." << std::endl;

	// Stop the running simulation
	projdyn_stop();

    // 将曲面网格转换为特征矩阵
    surface_mesh::Surface_mesh* mesh = viewer->getMesh();
    ProjDyn::Positions vertices(mesh->n_vertices(), 3);
    int j = 0;
    for (auto v : mesh->vertices()) {
        vertices.row(j) << (ProjDyn::Scalar1)mesh->position(v).x,
                (ProjDyn::Scalar1)mesh->position(v).y,
                (ProjDyn::Scalar1)mesh->position(v).z;
        ++j;
    }

    if (viewer->is_using_rods()) {
        MatrixXf* pos = viewer->getRodsPos();
        std::vector<Positions> rods;
        std::vector<size_t> rod_indices = viewer->getRodIndices();

        for (size_t ind = 0; ind < rod_indices.size(); ind++) {
            size_t rod_index = rod_indices.at(ind);
            size_t next_index = ind == rod_indices.size()-1 ? pos->cols() : rod_indices.at(ind+1);

            Positions rodPos;
            rodPos.resize(next_index - rod_index, 3);

            for (size_t i = rod_index; i < next_index; i++) {
                Vector3f p = pos->col(i);
                rodPos.row(i-rod_index) << p.x(), p.y(), p.z();
            }
            rods.push_back(rodPos);
        }

        sim.setRods(rods);
        sim.setMesh(vertices);
    } else {
        sim.setMesh(vertices);
    }

    return true;
}
void setup_demo_scene(Viewer* viewer, float ball_radius, size_t resolution, size_t num_rods, float seg_length, float hair_thickness) {
    viewer->addRodsOnBall(ball_radius, resolution, num_rods, seg_length, hair_thickness);
    default_constraints = false;
}

void init_projdyn_gui(Viewer* viewer) {
    //Simulation panel
	Window* pd_win = new Window(viewer, "Simulation");
   // pd_win->setPosition(Vector2i(15, 230));
	pd_win->setPosition(Vector2i(15, 15));
	pd_win->setLayout(new GroupLayout());

    PopupButton* setupPopup = new PopupButton(pd_win, "Create", ENTYPO_ICON_COG);
    Popup *popup = setupPopup->popup();
    popup->setLayout(new GroupLayout());

    new Label(popup, "Num of hairs");
    IntBox<int>* hairs_box = new IntBox<int>(popup, PROJDYN_NUM_HAIRS_INITIAL);
    hairs_box->setEditable(true);
    hairs_box->setCallback([](int num_hairs) {
        projdyn_num_hairs = num_hairs;
    });

     new Label(popup, "Num of Segments");
    IntBox<int>* res_box = new IntBox<int>(popup, PROJDYN_RES_INITIAL);
    res_box->setEditable(true);
    res_box->setCallback([](int res) {
        projdyn_res = res+1;
    });

   new Label(popup, "Segment Length");
    FloatBox<float>* seg_length_box = new FloatBox<float>(popup, PROJDYN_SEG_LENGTH_INITIAL);
    seg_length_box->setEditable(true);
    seg_length_box->setCallback([](float seg_length) {
        projdyn_seg_length = seg_length;
    });

    new Label(popup, "Hair Thickness");
    FloatBox<float>* hair_thickness_box = new FloatBox<float>(popup, PROJDYN_HAIR_THICKNESS_INITIAL);
    hair_thickness_box->setEditable(true);
    hair_thickness_box->setCallback([](float hair_thickness) {
        projdyn_hair_thickness = hair_thickness;
        });

    new Label(popup, "Head Radius");
    FloatBox<float>* radius_box = new FloatBox<float>(popup, PROJDYN_BALL_RADIUS_INITIAL);
    radius_box->setEditable(true);
    radius_box->setCallback([](float radius) {
        projdyn_ball_radius = radius;
    });

    Button* b = new Button(popup, "Set", ENTYPO_ICON_THUMBS_UP);
    b->setCallback([viewer,setupPopup]() {
        setup_demo_scene(viewer, projdyn_ball_radius, projdyn_res, projdyn_num_hairs, projdyn_seg_length, projdyn_hair_thickness);
        setupPopup->setPushed(false);
    });

    Button* runsim_b = new Button(pd_win, "Run Simulation", ENTYPO_ICON_CONTROLLER_PLAY);
	runsim_b->setCallback([viewer]() {
		projdyn_start(viewer);
	});

	Button* stopsim_b = new Button(pd_win, "Stop Simulation", ENTYPO_ICON_CONTROLLER_STOP);
	stopsim_b->setCallback([]() {
		projdyn_stop();
	});

	/*Button* reset_b = new Button(pd_win, "Reset Positions", ENTYPO_ICON_CW);
	reset_b->setCallback([viewer]() {
		bool was_active = projdyn_active;
		projdyn_stop();
		sim.resetPositions();
		if (was_active) {
			projdyn_start(viewer);
		} else {
			projdyn_upload_positions(viewer);
		}
	});
*/
    Button* wind_toggle_b = new Button(pd_win, "Toggle Wind");
    wind_toggle_b->setCallback([]() {
        sim.toggleWind();
    });

	Label* iterations_label = new Label(pd_win, "Num of Iterations: ");
	IntBox<int>* iterations_box = new IntBox<int>(pd_win, projdyn_num_iterations);
    iterations_box->setEditable(true);
    iterations_box->setCallback([viewer](int num_its) {
		projdyn_num_iterations = num_its;
	});

	viewer->performLayout();
}

bool projdyn_start(Viewer* viewer) {
	projdyn_stop();

	// Make sure the simulator is properly initialized
	if (!sim.isInitialized()) {
        auto t1 = std::chrono::high_resolution_clock::now();

        bool init = sim.initializeSystem(default_constraints);

        auto t2 = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

        std::cout << "System initialized in: " << duration << "s" << std::endl;

        if (!init)

            return false;

    }

	// 创建一个运行模拟的线程
    // 它调用一个函数，每 1000/PROJDYN_FPS 毫秒触发一次时间步长
	projdyn_active = true;
	projdyn_thread = std::thread(
		[viewer]() {
		std::chrono::milliseconds time(1000 / PROJDYN_FPS);
		while (projdyn_active) {
			std::this_thread::sleep_for(time);
			projdyn_update(viewer);
			glfwPostEmptyEvent();
		}
	}
	);

	num_iter = 0;

	return true;
}

void projdyn_stop() {
	if (projdyn_active) {
        projdyn_active = false;
        projdyn_thread.join();
    }
}

// 每帧都会被调用：
// 执行时间步并更新在着色器窗口中绘制的位置
bool projdyn_update(Viewer* viewer) {
	if (!sim.isInitialized()) return false;

	const size_t max_iter = 100;
	if (num_iter > max_iter) {
	    num_iter = 0;
        float avg = 0.0f;
        for (size_t i = 0; i < exec_time.size(); i++) {
            avg += exec_time.at(i);
        }
        avg /= exec_time.size();
        exec_time.clear();

        std::cout << "Average execution: " << avg / 1e6 << "s" << std::endl;
	}

	//测量执行时间
    auto t1 = std::chrono::high_resolution_clock::now();

    // 模拟一个时间步
    sim.step(projdyn_num_iterations);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    exec_time.push_back(duration);

    num_iter++;

    return projdyn_upload_positions(viewer);
}

// 提取位置，将它们转换为按列的浮点数三元组和
// 将它们上传到 OpenGL 缓冲区
bool projdyn_upload_positions(Viewer* viewer) {
	// In this function you need to extract the vertex positions of the simulation
	// and send them to the viewer.

	Positions* pos = sim.getPositions();

    // 这是使用3 x#Verts矩阵upload_pos完成的
    size_t num_verts = pos->rows();
    upload_pos.resize(3, num_verts);

    //复制矩阵的转置，因为存储顺序不同
    for (size_t i = 0; i < num_verts; i++) {
        upload_pos(0, i) = pos->coeff(i, 0);
        upload_pos(1, i) = pos->coeff(i, 1);
        upload_pos(2, i) = pos->coeff(i,2);
    }

    std::vector<ProjDyn::Index> rod_indices;
    Positions* rods_pos = sim.getRodsPositions();
    size_t num_pos = rods_pos->rows();
    upload_rods.resize(3, num_pos);

    //获得n*3的矩阵
    for (size_t i = 0; i < num_pos; i++) {
        upload_rods(0, i) = rods_pos->coeff(i, 0);
        upload_rods(1, i) = rods_pos->coeff(i, 1);
        upload_rods(2, i) = rods_pos->coeff(i,2);
    }

    upload_rods_tan.resize(3, num_pos);
    Positions* rods_tan = sim.getRodsTangents();

    for (size_t i = 0; i < num_pos; i++) {
        upload_rods_tan(0, i) = rods_tan->coeff(i, 0);
        upload_rods_tan(1, i) = rods_tan->coeff(i, 1);
        upload_rods_tan(2, i) = rods_tan->coeff(i,2);
    }

    upload_rods_norm.resize(3, num_pos);
    Positions* rods_normals = sim.getRodsNormals();

    for (size_t i = 0; i < num_pos; i++) {
        upload_rods_norm(0, i) = rods_normals->coeff(i, 0);
        upload_rods_norm(1, i) = rods_normals->coeff(i, 1);
        upload_rods_norm(2, i) = rods_normals->coeff(i,2);
    }

    rod_indices = sim.getRodIndices();

	// The matrix is sent to the viewer with this function
    viewer->updateShaderRods(upload_rods, upload_rods_tan, upload_rods_norm, rod_indices);
    viewer->updateShaderVertices(upload_pos);

	return true;
}

void projdyn_grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
	if (!sim.isInitialized() || !projdyn_active) return;

	sim.setGrab(grabbedVerts, grabPos);
}

void projdyn_release_grab() {
	if (!sim.isInitialized()) return;

	sim.releaseGrab();
}

#include <Novice.h>
#include <Matrix.h>

const char kWindowTitle[] = "GC2B_02_イトウ_シュンヤ";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Matrix matrix;

	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };
	Vector3 cameraPosition = { 0.0f,0.f,-6.0f };

	Matrix4x4 cameraMatrix = Matrix::MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.26f,0.0f,0.0f }, cameraPosition);
	Matrix4x4 viewMatrix = Matrix::Inverse(cameraMatrix);
	Matrix4x4 projectionMatrix = Matrix::MakePerspectiveFovMatrix(0.45f, float(1280) / float(720), 0.1f, 100.0f);
	Matrix4x4 viewportMatrix = Matrix::MakeViewportMatrix(0, 0, float(1280), float(720), 0.0f, 1.0f);
	Matrix4x4 viewProjectionMatrix = Matrix::Multiply(viewMatrix, projectionMatrix);

	//正射影ベクトルと最近接点
	Matrix::Segment segment{ 
		.origin{-0.7f, 0.3f, 0.0f},
		.diff  {2.0f, -0.5f, 0.0f}
	};
	//線分の描画
	Vector3 start = Matrix::Transform(Matrix::Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
	Vector3 end = Matrix::Transform(Matrix::Transform(Matrix::Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);
	uint32_t segmentColor = 0xFFFFFFFF;


	Matrix::AABB aabb{
		.min{-0.5f,-0.5f,-0.5f},
		.max{0.5f ,0.5f ,0.5f},
	};

	uint32_t aabbColor = 0xFFFFFFFF;

	bool iscollision = Matrix::IsCollision(aabb, segment);

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		cameraMatrix = Matrix::MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraPosition);
		viewMatrix = Matrix::Inverse(cameraMatrix);
		projectionMatrix = Matrix::MakePerspectiveFovMatrix(0.45f, float(1280) / float(720), 0.1f, 100.0f);
		viewportMatrix = Matrix::MakeViewportMatrix(0, 0, float(1280), float(720), 0.0f, 1.0f);
		viewProjectionMatrix = Matrix::Multiply(viewMatrix, projectionMatrix);

		start = Matrix::Transform(Matrix::Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
		end = Matrix::Transform(Matrix::Transform(Matrix::Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);

		iscollision = Matrix::IsCollision(aabb, segment);
		if (iscollision == true) {
			aabbColor = RED;
		}
		else { aabbColor = WHITE; }

		ImGui::Begin("Window");

		ImGui::DragFloat3("CameraTranslate", &cameraPosition.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);

		ImGui::DragFloat3("aabb.min", &aabb.min.x, 0.01f);
		ImGui::DragFloat3("aabb.max", &aabb.max.x, 0.01f);

		ImGui::DragFloat3("Segment origin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("Segment diff", &segment.diff.x, 0.01f);


		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		matrix.DrawGrid(viewProjectionMatrix, viewportMatrix);
		matrix.DrawAABB(aabb, viewProjectionMatrix, viewportMatrix, aabbColor);
		Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), segmentColor);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
